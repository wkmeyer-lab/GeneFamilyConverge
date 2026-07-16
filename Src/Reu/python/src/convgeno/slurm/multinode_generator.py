"""SLURM script generation for multi-node OrthoFinder execution.

The multi-node mode splits OrthoFinder into three chained jobs:

1. **Prepare** — runs ``orthofinder -op`` on a single node to format inputs
   and emit the list of DIAMOND/BLAST search commands.
2. **Search array** — a SLURM array job that distributes those search
   commands across multiple nodes.
3. **Resume** — runs ``orthofinder -b <workdir>`` on a single node once
   all search commands are done, performing clustering, MSA, and tree
   inference.
"""

from __future__ import annotations

import heapq
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from convgeno.slurm.config import (
    MultinodeConfig,
    PipelineConfig,
    normalize_optional_account,
)
from convgeno.slurm.discovery import (
    detect_max_array_size,
    detect_node_cpus,
    detect_qos_max_jobs,
)
from convgeno.slurm.runtime import CondaRuntimeConfig, render_conda_bootstrap

# POSIX ERE alternation matching a real OrthoFinder search command.
#
# Requirements baked in:
#   * The program name must be followed by a flag token (starting with
#     "-"). OrthoFinder's prose header lines like "diamond commands that
#     must be run" or "blastp commands that must be run" do NOT start
#     with a flag, so this filter rejects them. The previous loose
#     pattern ``^(diamond|blastp|makeblastdb)`` matched those headers
#     and the search array crashed trying to ``eval`` them.
#   * For ``diamond`` we additionally require the subcommand keyword
#     (``blastp`` or ``makedb``).
_COMMAND_START_ALT = (
    r"diamond[[:space:]]+(blastp|makedb)[[:space:]]+-"
    r"|blastp[[:space:]]+-"
    r"|makeblastdb[[:space:]]+-"
)

# Validation regex: matches a *cleaned* command line in $COMMANDS_FILE
# (leading whitespace already stripped by sed).
COMMAND_LINE_REGEX = f"^({_COMMAND_START_ALT})"

# Extraction regex: matches a raw log line, which OrthoFinder sometimes
# indents. Leading whitespace is stripped by the sed stage of the
# extraction pipeline before the line is written to $COMMANDS_FILE.
COMMAND_LINE_REGEX_EXTRACT = f"^[[:space:]]*({_COMMAND_START_ALT})"


def _build_sbatch_header(config: PipelineConfig, overrides: dict) -> str:
    """Build SBATCH header lines from ``config.slurm`` with field overrides.

    Returns the lines joined into a single ``\\n``-separated string.
    """
    lines: dict[str, str | None] = {
        "--partition": config.slurm.partition,
        "--nodes": str(config.slurm.nodes),
        "--ntasks": str(config.slurm.ntasks),
        "--cpus-per-task": str(config.slurm.cpus_per_task),
        "--time": config.slurm.time_limit,
        "--mem": config.slurm.mem,
        "--mem-per-cpu": config.slurm.mem_per_cpu,
        "--output": config.slurm.output_pattern,
        "--error": config.slurm.error_pattern,
        "--export": "ALL",
    }
    account = normalize_optional_account(config.slurm.account)
    if account is not None:
        lines["--account"] = account
    if config.slurm.mail_user is not None:
        lines["--mail-user"] = config.slurm.mail_user
        lines["--mail-type"] = config.slurm.mail_type

    lines.update(overrides)

    formatted = [
        f"#SBATCH {flag}={value}"
        for flag, value in lines.items()
        if value is not None
    ]
    for arg in config.slurm.extra_sbatch_args:
        formatted.append(f"#SBATCH {arg}")
    return "\n".join(formatted)


def _validate_has_orthofinder(config: PipelineConfig) -> None:
    if config.orthofinder is None:
        raise ValueError(
            "OrthoFinder settings not found in pipeline config. "
            "Add an 'orthofinder' section to pipeline_config.yaml."
        )


def _resolve_runtime(
    config: PipelineConfig,
    runtime: CondaRuntimeConfig | None,
) -> CondaRuntimeConfig:
    """Resolve runtime from explicit parameter or config, raising if absent."""
    resolved = runtime if runtime is not None else config.runtime
    if resolved is None:
        raise ValueError(
            "No runtime configuration provided. "
            "Run 'convgeno init' to detect and store conda paths."
        )
    return resolved


# ============================================================
# Prepare-phase resource derivation
# ============================================================
# Prepare parses every proteome and builds one search database per species,
# so its cost scales with the TOTAL proteome volume, and its peak memory with
# the LARGEST single proteome (diamond makedb loads one proteome at a time).
# These are conservative, generalizable heuristics computed at generate time
# from the input FASTA sizes -- not values fitted to any particular run.
_PREPARE_CPU_CAP = 16
_PREPARE_MEM_BASE_MB = 2048  # interpreter + OrthoFinder overhead
_PREPARE_MEM_PER_TOTAL_MB = 2  # x total proteome MB (all-gene ID parsing)
_PREPARE_MEM_PER_LARGEST_MB = 4  # x largest proteome MB (makedb peak)
_PREPARE_MEM_FLOOR_MB = 4096
_PREPARE_MEM_ROUND_MB = 1024
_PREPARE_TIME_BASE_SEC = 600  # fixed startup / I/O overhead
_PREPARE_TIME_PER_TOTAL_MB_SEC = 3
_PREPARE_TIME_MARGIN = 1.5
_PREPARE_TIME_MIN_SEC = 1800  # 30 min floor
_PREPARE_TIME_MAX_SEC = 14400  # 4 h cap
_PREPARE_FALLBACK_MEM_MB = _PREPARE_MEM_FLOOR_MB
_PREPARE_FALLBACK_TIME = "02:00:00"

_FASTA_EXTENSIONS = (".fa", ".fasta", ".faa")


def _seconds_to_hms(seconds: int) -> str:
    """Format whole seconds as SLURM ``HH:MM:SS``."""
    return f"{seconds // 3600:02d}:{(seconds % 3600) // 60:02d}:{seconds % 60:02d}"


def derive_prepare_cpus(cpus_per_task: int) -> int:
    """Prepare CPU count = ``min(16, cores)``.

    ``cpus_per_task`` is the configured per-node allocation, which ``init``
    sets to (partition physical cores - 4). Prepare is a light parse+makedb
    job, so it is capped at 16 cores (and uses fewer on small nodes) -- enough
    for makedb without reserving a whole fat node, which keeps it easy to
    backfill.
    """
    return min(_PREPARE_CPU_CAP, cpus_per_task)


def derive_prepare_memory_mb(total_bytes: int, largest_bytes: int) -> int:
    """Derive prepare memory (MB) from proteome volume.

    Peak memory tracks the largest single proteome (makedb) plus a term for
    OrthoFinder's all-gene ID parsing (total volume), over a base overhead;
    floored and rounded up.
    """
    total_mb = total_bytes / (1024 * 1024)
    largest_mb = largest_bytes / (1024 * 1024)
    raw = (
        _PREPARE_MEM_BASE_MB
        + _PREPARE_MEM_PER_TOTAL_MB * total_mb
        + _PREPARE_MEM_PER_LARGEST_MB * largest_mb
    )
    rounded = math.ceil(raw / _PREPARE_MEM_ROUND_MB) * _PREPARE_MEM_ROUND_MB
    return max(_PREPARE_MEM_FLOOR_MB, int(rounded))


def derive_prepare_walltime(total_bytes: int) -> str:
    """Derive prepare walltime (HH:MM:SS) from total proteome volume.

    Parsing all proteomes and building the databases scales with total
    sequence volume; a safety margin is applied and the result clamped to a
    sane [30 min, 4 h] range.
    """
    total_mb = total_bytes / (1024 * 1024)
    seconds = _PREPARE_TIME_BASE_SEC + _PREPARE_TIME_PER_TOTAL_MB_SEC * total_mb
    seconds = int(seconds * _PREPARE_TIME_MARGIN)
    seconds = max(_PREPARE_TIME_MIN_SEC, min(_PREPARE_TIME_MAX_SEC, seconds))
    return _seconds_to_hms(seconds)


def _proteome_volume(input_dir: str) -> tuple[int, int, int]:
    """Return ``(total_bytes, largest_bytes, n_files)`` for FASTA proteomes.

    Returns ``(0, 0, 0)`` when the directory is absent or unreadable (e.g.
    generating scripts off-cluster before the data is staged), so callers fall
    back to safe defaults.
    """
    try:
        directory = Path(input_dir)
        if not directory.is_dir():
            return (0, 0, 0)
        sizes = [
            entry.stat().st_size
            for entry in directory.iterdir()
            if entry.is_file() and entry.suffix.lower() in _FASTA_EXTENSIONS
        ]
    except OSError:
        return (0, 0, 0)
    if not sizes:
        return (0, 0, 0)
    return (sum(sizes), max(sizes), len(sizes))


def generate_prepare_script(
    config: PipelineConfig,
    runtime: CondaRuntimeConfig | None = None,
    sizing: SearchArraySizing | None = None,
) -> str:
    """Generate SLURM script for OrthoFinder's prepare phase (``-op``).

    Runs on a single node. Captures the DIAMOND/BLAST commands to a file, then
    builds the balanced per-task search manifest so the array job can consume
    it. The bucket count ``T`` is taken from ``sizing`` (derived once and shared
    with the search-array script so the manifest and ``--array`` width agree).
    """
    _validate_has_orthofinder(config)
    assert config.orthofinder is not None  # for type checker
    resolved_runtime = _resolve_runtime(config, runtime)
    if sizing is None:
        sizing = derive_search_array_sizing(config)

    # Derive prepare resources from the proteome volume at generate time
    # (SBATCH headers are static). Fall back to safe defaults when the inputs
    # are not readable yet (e.g. generating off-cluster before data is staged).
    prepare_cpus = derive_prepare_cpus(config.slurm.cpus_per_task)
    total_bytes, largest_bytes, _ = _proteome_volume(config.orthofinder.input_dir)
    if total_bytes > 0:
        prepare_mem = f"{derive_prepare_memory_mb(total_bytes, largest_bytes)}M"
        prepare_time = derive_prepare_walltime(total_bytes)
    else:
        prepare_mem = f"{_PREPARE_FALLBACK_MEM_MB}M"
        prepare_time = _PREPARE_FALLBACK_TIME

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header = _build_sbatch_header(
        config,
        overrides={
            "--job-name": "convgeno_of_prepare",
            "--nodes": "1",
            "--ntasks": "1",
            "--cpus-per-task": str(prepare_cpus),
            "--time": prepare_time,
            "--mem": prepare_mem,
            "--mem-per-cpu": None,
        },
    )
    bootstrap_block = render_conda_bootstrap(resolved_runtime)

    script = f"""\
#!/bin/bash
# ============================================================
# OrthoFinder prepare phase (-op) — multi-node mode
# Generated by convgeno on {timestamp}
# ============================================================

{header}

set -euo pipefail

{bootstrap_block}

echo "============================================================"
echo "convgeno: OrthoFinder prepare phase started"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $HOSTNAME"
echo "Start:     $(date)"
echo "============================================================"

START_SECONDS=$SECONDS

INPUT_DIR="{config.orthofinder.input_dir}"
OUTPUT_DIR="{config.orthofinder.output_dir}"
SEARCH_PROGRAM="{config.orthofinder.sequence_search}"

# All auxiliary files (log, commands list, WorkingDirectory pointer)
# live in the PARENT directory, not in $OUTPUT_DIR itself. OrthoFinder
# refuses to run if its -o output directory already exists, so we must
# not pre-create $OUTPUT_DIR and must not redirect any output into it
# before OrthoFinder has created it.
OUTPUT_PARENT="$(dirname "$OUTPUT_DIR")"
RUN_NAME="$(basename "$OUTPUT_DIR")"

mkdir -p "$OUTPUT_PARENT"

if [ -e "$OUTPUT_DIR" ]; then
    echo "ERROR: OrthoFinder output directory already exists: $OUTPUT_DIR"
    echo "Choose a fresh output_dir in the config."
    exit 1
fi

PREPARE_LOG="$OUTPUT_PARENT/${{RUN_NAME}}_prepare_full_stdout.log"
COMMANDS_FILE="$OUTPUT_PARENT/${{RUN_NAME}}_diamond_commands.txt"
WORK_DIR_FILE="$OUTPUT_PARENT/${{RUN_NAME}}_working_dir_path.txt"

# Run OrthoFinder prepare phase. -op stops after writing the search
# commands and exits without running them. stdout is redirected to a
# log file in $OUTPUT_PARENT (which already exists), NOT in $OUTPUT_DIR.
orthofinder -f "$INPUT_DIR" -o "$OUTPUT_DIR" -op -S "$SEARCH_PROGRAM" > "$PREPARE_LOG" 2>&1

PREPARE_EXIT=$?
if [ $PREPARE_EXIT -ne 0 ]; then
    echo "ERROR: OrthoFinder prepare phase failed with exit code $PREPARE_EXIT"
    cat "$PREPARE_LOG"
    exit $PREPARE_EXIT
fi

# Extract search commands from OrthoFinder stdout. The regex requires a
# subcommand keyword (blastp/makedb) followed by a flag token (-...), so
# prose header lines like "diamond commands that must be run" are
# excluded. Leading whitespace is then stripped because OrthoFinder
# sometimes indents commands.
grep -E '{COMMAND_LINE_REGEX_EXTRACT}' "$PREPARE_LOG" \\
    | sed 's/^[[:space:]]*//' \\
    > "$COMMANDS_FILE" || true

# Validation: bail out *now* if the command file is empty or contains a
# line that the search-array job cannot safely `eval`. The previous bug
# (a header line in the file) escaped detection because the file was
# non-empty and counted line-by-line — search task 0 then tried to run
# `diamond commands that must be run` and crashed DIAMOND.
if [ ! -s "$COMMANDS_FILE" ]; then
    echo "ERROR: No valid DIAMOND/search commands were extracted."
    echo "Check prepare log: $PREPARE_LOG"
    exit 1
fi

BAD_LINES=$(grep -nEv '{COMMAND_LINE_REGEX}' "$COMMANDS_FILE" || true)

if [ -n "$BAD_LINES" ]; then
    echo "ERROR: Invalid lines found in command file:"
    echo "$BAD_LINES"
    echo "Command file: $COMMANDS_FILE"
    echo "Prepare log: $PREPARE_LOG"
    exit 1
fi

NUM_COMMANDS=$(wc -l < "$COMMANDS_FILE")
echo "Extracted $NUM_COMMANDS search commands to $COMMANDS_FILE"

# Locate the WorkingDirectory OrthoFinder created and persist its path
# for the resume script.
WORK_DIR=$(find "$OUTPUT_DIR" -maxdepth 2 -name "WorkingDirectory" -type d | head -1)

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: Could not find OrthoFinder WorkingDirectory under $OUTPUT_DIR"
    echo "Check prepare log: $PREPARE_LOG"
    exit 1
fi

echo "$WORK_DIR" > "$WORK_DIR_FILE"

# ------------------------------------------------------------
# Detect the OrthoFinder database-build behaviour, split the emitted
# commands, and verify completeness (databases = n, searches = n^2).
#
# Some OrthoFinder builds create the DIAMOND/BLAST databases themselves
# during -op and emit only the search (blastp) commands; others emit the
# DB-build (makedb/makeblastdb) commands for us to run. The -op run above
# IS the probe (it never runs the searches), so we inspect what it produced
# -- no separate throwaway run is needed.
# ------------------------------------------------------------
DB_MODE_FILE="$OUTPUT_PARENT/${{RUN_NAME}}_db_mode.txt"
DB_BUILD_COMMANDS_FILE="$OUTPUT_PARENT/${{RUN_NAME}}_db_build_commands.txt"
SEARCH_COMMANDS_FILE="$OUTPUT_PARENT/${{RUN_NAME}}_search_commands.txt"
SPECIES_IDS_FILE="$WORK_DIR/SpeciesIDs.txt"

# Ground truth for completeness: one active line per species in
# SpeciesIDs.txt (removed species are commented with '#' and excluded).
if [ ! -f "$SPECIES_IDS_FILE" ]; then
    echo "ERROR: SpeciesIDs.txt not found: $SPECIES_IDS_FILE" >&2
    exit 1
fi
N_SPECIES=$(grep -cE '^[0-9]+:' "$SPECIES_IDS_FILE" || true)
if [ "$N_SPECIES" -lt 1 ]; then
    echo "ERROR: Could not read a species count from $SPECIES_IDS_FILE" >&2
    exit 1
fi
EXPECTED_SEARCHES=$(( N_SPECIES * N_SPECIES ))

# Split the extracted commands: DB-build (makedb/makeblastdb) vs search
# (blastp). The search array consumes SEARCH_COMMANDS_FILE (blastp only), so
# it never re-runs makedb and its length is exactly n^2.
grep -E 'makedb|makeblastdb' "$COMMANDS_FILE" > "$DB_BUILD_COMMANDS_FILE" || true
grep -Ev 'makedb|makeblastdb' "$COMMANDS_FILE" > "$SEARCH_COMMANDS_FILE" || true
DB_BUILD_COUNT=$(wc -l < "$DB_BUILD_COMMANDS_FILE")
SEARCH_COUNT=$(wc -l < "$SEARCH_COMMANDS_FILE")
# maxdepth 1 only: OrthoFinder's startup self-test leaves a throwaway test DB
# under $WORK_DIR/dependencies/, which must never be counted as a species DB.
DB_FILE_COUNT=$(find "$WORK_DIR" -maxdepth 1 \\( -name '*.dmnd' -o -name '*.phr' -o -name '*.pin' -o -name '*.psq' \\) 2>/dev/null | wc -l)

echo ""
echo "---- OrthoFinder DB behaviour + completeness ----"
echo "Species (n):               $N_SPECIES"
echo "DB-build commands emitted: $DB_BUILD_COUNT"
echo "Prebuilt DB files present: $DB_FILE_COUNT"
echo "Search (blastp) commands:  $SEARCH_COUNT  (expected n^2 = $EXPECTED_SEARCHES)"

# Determine the DB mode from two independent signals.
if [ "$DB_BUILD_COUNT" -gt 0 ]; then
    OF_DB_MODE="emit_build_commands"
    echo "Detected DB mode: EMIT_BUILD_COMMANDS (this build does not create the DBs itself)."
elif [ "$DB_FILE_COUNT" -gt 0 ]; then
    OF_DB_MODE="self_built"
    echo "Detected DB mode: SELF_BUILT (-op already created the databases)."
else
    echo "ERROR: Could not determine how OrthoFinder handled the search databases." >&2
    echo "  No makedb/makeblastdb build commands were emitted, and no" >&2
    echo "  *.dmnd/*.phr/*.pin/*.psq database files were found under:" >&2
    echo "  $WORK_DIR" >&2
    echo "  Prepare log: $PREPARE_LOG" >&2
    exit 1
fi
echo "$OF_DB_MODE" > "$DB_MODE_FILE"

# (1) Searches must be exactly n^2 (all ordered species pairs, incl. self).
if [ "$SEARCH_COUNT" -ne "$EXPECTED_SEARCHES" ]; then
    echo "ERROR: Incomplete search list -- found $SEARCH_COUNT blastp commands," >&2
    echo "  expected n^2 = $EXPECTED_SEARCHES for $N_SPECIES species." >&2
    echo "  Proceeding would silently corrupt orthogroups; aborting. See $PREPARE_LOG." >&2
    exit 1
fi

# (2) Databases: one per species (n). Build them here when OrthoFinder didn't.
if [ "$OF_DB_MODE" = "emit_build_commands" ]; then
    if [ "$DB_BUILD_COUNT" -ne "$N_SPECIES" ]; then
        echo "ERROR: Expected $N_SPECIES database-build commands (one per species)," >&2
        echo "  but found $DB_BUILD_COUNT. Aborting. See $PREPARE_LOG." >&2
        exit 1
    fi
    echo "Building $DB_BUILD_COUNT search databases (one per species) on this node..."
    DB_BUILD_FAILED=0
    while IFS= read -r BUILD_CMD; do
        if [ -z "$BUILD_CMD" ]; then
            continue
        fi
        echo "[makedb] $BUILD_CMD"
        if ! eval "$BUILD_CMD"; then
            echo "WARNING: database-build command failed: $BUILD_CMD" >&2
            DB_BUILD_FAILED=$(( DB_BUILD_FAILED + 1 ))
        fi
    done < "$DB_BUILD_COMMANDS_FILE"
    if [ "$DB_BUILD_FAILED" -gt 0 ]; then
        echo "ERROR: $DB_BUILD_FAILED database-build command(s) failed; aborting." >&2
        exit 1
    fi
fi

# (3) Verify the databases now exist -- one per species, at the TOP LEVEL of
# the WorkingDirectory. Check each species' DB by name (maxdepth 1): the real
# per-species DBs live directly in the WorkingDirectory, while OrthoFinder's
# startup self-test leaves a throwaway test DB under $WORK_DIR/dependencies/
# that must be ignored.
if [ "$SEARCH_PROGRAM" = "diamond" ]; then
    MISSING_DBS=0
    for SPECIES_ID in $(grep -oE '^[0-9]+' "$SPECIES_IDS_FILE" || true); do
        if [ ! -f "$WORK_DIR/diamondDBSpecies${{SPECIES_ID}}.dmnd" ]; then
            echo "ERROR: Missing DIAMOND database: diamondDBSpecies${{SPECIES_ID}}.dmnd" >&2
            MISSING_DBS=$(( MISSING_DBS + 1 ))
        fi
    done
    if [ "$MISSING_DBS" -gt 0 ]; then
        echo "ERROR: $MISSING_DBS of $N_SPECIES DIAMOND databases missing" >&2
        echo "  (top-level diamondDBSpecies<id>.dmnd in $WORK_DIR). Aborting. See $PREPARE_LOG." >&2
        exit 1
    fi
    echo "Verified $N_SPECIES DIAMOND databases (per-species, top-level)."
else
    OTHER_DB_COUNT=$(find "$WORK_DIR" -maxdepth 1 \\( -name '*.phr' -o -name '*.pin' -o -name '*.psq' -o -name '*.pdb' \\) 2>/dev/null | wc -l)
    if [ "$OTHER_DB_COUNT" -lt 1 ]; then
        echo "ERROR: No search databases found (top-level) in $WORK_DIR for '$SEARCH_PROGRAM'." >&2
        exit 1
    fi
    echo "Verified search databases present ($OTHER_DB_COUNT files) for '$SEARCH_PROGRAM'."
fi

echo "Search commands (blastp only, n^2=$EXPECTED_SEARCHES) -> $SEARCH_COMMANDS_FILE"
echo "-------------------------------------------------"

# ------------------------------------------------------------
# Build the balanced per-task search manifest (LPT) for the array job.
# T (task count) is fixed at generate time so the array width matches exactly.
# ------------------------------------------------------------
MANIFEST_DIR="$OUTPUT_PARENT/${{RUN_NAME}}{_SEARCH_MANIFEST_DIRNAME_SUFFIX}"
SEARCH_TASKS={sizing.tasks}
echo "Building balanced search manifest: $SEARCH_TASKS tasks -> $MANIFEST_DIR"
if ! python -m convgeno.slurm.build_search_manifest \\
        --search-commands "$SEARCH_COMMANDS_FILE" \\
        --work-dir "$WORK_DIR" \\
        --manifest-dir "$MANIFEST_DIR" \\
        --tasks "$SEARCH_TASKS"; then
    echo "ERROR: failed to build the search manifest; aborting prepare." >&2
    exit 1
fi
echo "Search manifest ready in $MANIFEST_DIR"

ELAPSED=$(( SECONDS - START_SECONDS ))
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "============================================================"
echo "convgeno: prepare phase finished"
echo "End:       $(date)"
echo "Duration:  ${{HOURS}}h ${{MINUTES}}m ${{SECS}}s"
echo "============================================================"
exit 0
"""
    return script


# ============================================================
# Search-array load balancing (LPT bucketing)
# ============================================================
# The prepare phase emits exactly n^2 `diamond blastp` commands (one per
# ordered species pair, self-searches included). The search array must
# distribute them across T tasks so every task finishes at about the same
# time -- a multiprocessor makespan-minimisation problem.
#
# We use greedy Longest-Processing-Time (LPT): weight each command by its
# estimated cost, sort descending, then assign each to the currently
# least-loaded bucket. Graham's bound guarantees makespan <=
# (4/3 - 1/(3T)) * OPT, so buckets are provably near-balanced -- unlike v1's
# cost-blind contiguous chunks, where one all-big-searches chunk set the whole
# array's wall time.
#
# Cost model (a proxy for diamond blastp runtime): cost(i, j) = |S_i| * |S_j|,
# the product of the query and target proteome sizes in bytes. Both drive the
# alignment work, so self-searches and big x big pairs dominate.
#
# All arithmetic is integer on purpose: |S_i| * |S_j| for a large proteome
# exceeds float's exact-integer range (2**53), so floats would lose precision
# and make bucket sums -- and therefore tie-breaks -- non-deterministic.


@dataclass(frozen=True)
class SearchCommand:
    """One emitted ``diamond blastp`` command, with its LPT cost inputs.

    ``line_index`` is the command's 0-based position in the search-commands
    file (``_search_commands.txt``); a per-task manifest is just a list of
    these indices. ``query_bytes`` / ``db_bytes`` are the query (``S_i``) and
    target-database (``S_j``) proteome sizes in bytes. ``db_bytes`` is tracked
    separately because it -- not the query -- drives a diamond run's peak
    memory: diamond loads the target database into RAM and streams the query.
    """

    line_index: int
    query_bytes: int
    db_bytes: int

    @property
    def cost(self) -> int:
        """LPT weight ``|S_i| * |S_j|`` (bytes^2), a proxy for runtime."""
        return self.query_bytes * self.db_bytes


@dataclass(frozen=True)
class LPTResult:
    """Result of :func:`lpt_partition`.

    ``buckets[b]`` holds the ``line_index`` values assigned to task ``b``
    (ascending); ``bucket_costs[b]`` is that bucket's summed cost (parallel to
    ``buckets``). ``max_bucket_cost`` -- the makespan-determining bucket --
    later sizes the array's per-task ``--time``; ``mem_determinant_bytes`` --
    the largest target database over *all* commands -- later sizes its per-task
    ``--mem``. Neither is converted to SLURM units here: that needs C/p, base,
    a throughput constant, and margins from the later sizing sub-steps.
    """

    buckets: list[list[int]]
    bucket_costs: list[int]
    max_bucket_cost: int
    mem_determinant_bytes: int


def lpt_partition(
    commands: list[SearchCommand],
    num_buckets: int,
) -> LPTResult:
    """Distribute ``commands`` into ``num_buckets`` cost-balanced buckets (LPT).

    Sorts commands by cost descending (ties broken by ``line_index`` for
    determinism) and greedily assigns each to the least-loaded bucket via a
    min-heap. The first ``num_buckets`` commands therefore seed distinct empty
    buckets, so no bucket is empty when ``num_buckets <= len(commands)``.

    ``mem_determinant_bytes`` is the largest ``db_bytes`` over *all* commands,
    not just the biggest bucket: a SLURM array applies one ``--mem`` to every
    task, so it must cover whichever task ends up holding the largest database.

    Raises ``ValueError`` if ``commands`` is empty or ``num_buckets < 1``.
    """
    if num_buckets < 1:
        raise ValueError(f"num_buckets must be >= 1, got {num_buckets}")
    if not commands:
        raise ValueError("commands must be non-empty")

    # Sort by cost descending; line_index breaks ties for a stable order.
    ordered = sorted(commands, key=lambda c: (-c.cost, c.line_index))

    buckets: list[list[int]] = [[] for _ in range(num_buckets)]
    bucket_costs = [0] * num_buckets
    # Min-heap of (current_total_cost, bucket_index). On ties the lowest
    # bucket_index pops first, keeping the assignment fully deterministic.
    heap = [(0, b) for b in range(num_buckets)]
    heapq.heapify(heap)

    for cmd in ordered:
        _, b = heapq.heappop(heap)
        buckets[b].append(cmd.line_index)
        bucket_costs[b] += cmd.cost
        heapq.heappush(heap, (bucket_costs[b], b))

    for bucket in buckets:
        bucket.sort()

    return LPTResult(
        buckets=buckets,
        bucket_costs=bucket_costs,
        max_bucket_cost=max(bucket_costs),
        mem_determinant_bytes=max(cmd.db_bytes for cmd in commands),
    )


# ------------------------------------------------------------
# Building the SearchCommand list from real prepare output
# ------------------------------------------------------------
# `lpt_partition` consumes SearchCommands; these helpers construct them from
# the two artefacts the prepare job leaves in the WorkingDirectory: the
# blastp-only `_search_commands.txt` (one command per line) and the per-species
# `Species{id}.fa` proteomes (their byte sizes are the cost/memory proxy).
#
# The species pair is read from the `Blast{i}_{j}` output token rather than by
# matching an exact flag layout: OrthoFinder's own `-b` resume looks for
# `Blast{i}_{j}.txt.gz`, so that naming is a hard, version-stable contract,
# whereas the surrounding diamond flags vary across versions. `i` is the query
# species, `j` is the target-database species (cost = |S_i|*|S_j|, memory ~
# |S_j|). The query (`Species{i}.fa`) and database (`diamondDBSpecies{j}`)
# tokens are cross-checked when present, and parsing fails loud otherwise.

# Match `Blast<i>_<j>` anywhere on the line (the -o output path).
_BLAST_PAIR_RE = re.compile(r"Blast(\d+)_(\d+)")
# Match a standalone `Species<i>.fa` query token; the lookbehind stops it from
# matching the `Species` inside `diamondDBSpecies<j>`.
_QUERY_SPECIES_RE = re.compile(r"(?<![A-Za-z])Species(\d+)\.fa")
# Match the diamond database token `diamondDBSpecies<j>` (with or without the
# `.dmnd` extension).
_DB_SPECIES_RE = re.compile(r"diamondDBSpecies(\d+)")
# Match an exact `Species<id>.fa` filename (for reading proteome sizes).
_SPECIES_FASTA_FILE_RE = re.compile(r"^Species(\d+)\.fa$")


def _parse_species_pair(command: str) -> tuple[int, int]:
    """Return ``(query_species_id, db_species_id)`` for one search command.

    Primary signal is the ``Blast{i}_{j}`` output token; the query
    (``Species{i}.fa``) and database (``diamondDBSpecies{j}``) tokens are used
    to cross-check it, or as a fallback when it is absent. Raises
    ``ValueError`` when the pair cannot be determined or the tokens disagree.
    """
    blast = _BLAST_PAIR_RE.search(command)
    query = _QUERY_SPECIES_RE.search(command)
    db = _DB_SPECIES_RE.search(command)

    if blast is not None:
        i, j = int(blast.group(1)), int(blast.group(2))
        if query is not None and int(query.group(1)) != i:
            raise ValueError(
                f"Query token Species{query.group(1)}.fa disagrees with output "
                f"Blast{i}_{j} in command: {command!r}"
            )
        if db is not None and int(db.group(1)) != j:
            raise ValueError(
                f"Database token diamondDBSpecies{db.group(1)} disagrees with "
                f"output Blast{i}_{j} in command: {command!r}"
            )
        return i, j

    # No Blast{i}_{j} token: fall back to the query + database tokens.
    if query is not None and db is not None:
        return int(query.group(1)), int(db.group(1))

    raise ValueError(
        "Could not determine the (query, database) species pair from search "
        f"command: {command!r}"
    )


def read_species_fasta_sizes(work_dir: Path) -> dict[int, int]:
    """Map species id -> ``Species{id}.fa`` byte size from a WorkingDirectory.

    OrthoFinder writes one ``Species{id}.fa`` per species during ``-op``; its
    byte size is the proteome-volume proxy the cost model uses. Ignores
    everything else (``SpeciesIDs.txt``, ``diamondDBSpecies*.dmnd``, ...).
    Raises ``FileNotFoundError`` when no such files exist.
    """
    sizes: dict[int, int] = {}
    for entry in Path(work_dir).iterdir():
        match = _SPECIES_FASTA_FILE_RE.match(entry.name)
        if match is not None and entry.is_file():
            sizes[int(match.group(1))] = entry.stat().st_size
    if not sizes:
        raise FileNotFoundError(
            f"No Species<id>.fa proteome files found in WorkingDirectory: "
            f"{work_dir}"
        )
    return sizes


def build_search_commands(
    command_lines: list[str],
    species_sizes: Mapping[int, int],
) -> list[SearchCommand]:
    """Build the LPT input list from emitted commands + proteome sizes.

    ``command_lines`` are the lines of ``_search_commands.txt`` (blastp only,
    one command each); ``line_index`` is each command's 0-based position, which
    the per-task manifest records and the array task uses to fetch the command.
    ``species_sizes`` maps species id -> ``Species{id}.fa`` byte size.

    Raises ``ValueError`` on a blank line, an unparseable command, or a species
    with no known size.
    """
    commands: list[SearchCommand] = []
    for line_index, raw in enumerate(command_lines):
        line = raw.strip()
        if not line:
            raise ValueError(
                f"Blank line at index {line_index} in the search-commands "
                "list; expected exactly one blastp command per line."
            )
        i, j = _parse_species_pair(line)
        for role, species in (("query", i), ("database", j)):
            if species not in species_sizes:
                raise ValueError(
                    f"No proteome size for {role} species {species} "
                    f"(command index {line_index}); known species: "
                    f"{sorted(species_sizes)}."
                )
        commands.append(
            SearchCommand(
                line_index=line_index,
                query_bytes=species_sizes[i],
                db_bytes=species_sizes[j],
            )
        )
    return commands


# ------------------------------------------------------------
# Per-task manifests (LPT buckets -> one file per array task)
# ------------------------------------------------------------
# lpt_partition returns buckets of line indices; the search array runs task
# `$SLURM_ARRAY_TASK_ID`, which reads its own manifest file and runs each line.
# Each manifest holds the resolved `diamond blastp` commands (copied from
# `_search_commands.txt`), so a task is self-contained: no cross-referencing
# back into the commands file, and the idempotent skip parses Blast{i}_{j}
# straight from the line it is about to run. `search_task_manifest_name` is the
# naming contract shared with the search-array script generator.

_SEARCH_MANIFEST_PREFIX = "search_task_"
_SEARCH_MANIFEST_SUFFIX = ".txt"
# Manifest directory basename suffix (appended to RUN_NAME). Shared contract
# between the prepare glue (writes here) and the search array (reads here).
_SEARCH_MANIFEST_DIRNAME_SUFFIX = "_search_manifests"


def search_task_manifest_name(task_id: int) -> str:
    """Manifest filename for array task ``task_id`` (== ``SLURM_ARRAY_TASK_ID``)."""
    return f"{_SEARCH_MANIFEST_PREFIX}{task_id}{_SEARCH_MANIFEST_SUFFIX}"


def write_task_manifests(
    buckets: list[list[int]],
    command_lines: list[str],
    manifest_dir: Path,
) -> list[Path]:
    """Write one manifest file per LPT bucket into ``manifest_dir``.

    ``buckets[t]`` holds line indices into ``command_lines`` (the lines of
    ``_search_commands.txt`` that produced the LPT partition); task ``t``'s file
    (``search_task_{t}.txt``) receives those commands, one per line. A file is
    written for **every** task id ``0..T-1`` -- empty buckets produce empty
    files -- so no array task ever reads a missing manifest. Returns the written
    paths in task order.

    Raises ``ValueError`` if a bucket references a line index outside
    ``command_lines`` (a caller contract violation: ``command_lines`` must be
    the same list ``build_search_commands`` indexed).
    """
    n = len(command_lines)
    manifest_dir = Path(manifest_dir)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []
    for task_id, bucket in enumerate(buckets):
        resolved: list[str] = []
        for idx in bucket:
            if not 0 <= idx < n:
                raise ValueError(
                    f"Task {task_id} references line index {idx}, out of range "
                    f"for {n} commands."
                )
            resolved.append(command_lines[idx].strip())
        path = manifest_dir / search_task_manifest_name(task_id)
        body = "\n".join(resolved)
        path.write_text(body + "\n" if body else "", encoding="utf-8")
        paths.append(path)
    return paths


# ------------------------------------------------------------
# Search-array sizing knobs (C, W, K, T)
# ------------------------------------------------------------
# C = cores per array task. Taking ~1/3 of the BIGGEST node makes each task a
# backfillable fraction rather than a whole-node reservation; capping at
# (SMALLEST node - 1) keeps a task schedulable on ANY node in the partition,
# maximising backfill. Both node counts come from detect_node_cpus()
# (discovery.py), which reports logical CPUs -- the unit SLURM --cpus-per-task
# uses. C is user-overridable (later, via MultinodeConfig.search_cpus).

_SEARCH_CPU_NODE_FRACTION = 3  # C targets ~1/3 of the biggest node
_SEARCH_CPU_SMALLEST_HEADROOM = 1  # leave >=1 core free on the smallest node


def derive_search_cpus(
    max_cpus_per_node: int,
    min_cpus_per_node: int,
    override: int | None = None,
) -> int:
    """Cores per search-array task (C), clamped to >= 1.

    ``C = min(max_cpus_per_node // 3, min_cpus_per_node - 1)``. The first term
    keeps a task to a backfillable fraction of the largest node; the second
    guarantees it still fits (with a core to spare) on the smallest node, so it
    can land anywhere in the partition. A positive ``override`` (from
    ``MultinodeConfig.search_cpus``) wins outright. When cluster discovery is
    unavailable (both counts 0), returns 1 -- regenerate on-cluster for a
    meaningful value.
    """
    if override is not None:
        if override < 1:
            raise ValueError(f"search_cpus override must be >= 1, got {override}")
        return override
    third_of_biggest = max_cpus_per_node // _SEARCH_CPU_NODE_FRACTION
    fits_smallest = min_cpus_per_node - _SEARCH_CPU_SMALLEST_HEADROOM
    return max(1, min(third_of_biggest, fits_smallest))


def derive_search_cpus_from_discovery(
    node_cpus: Mapping[str, int],
    override: int | None = None,
) -> int:
    """``derive_search_cpus`` fed from a ``detect_node_cpus`` result dict."""
    return derive_search_cpus(
        max_cpus_per_node=node_cpus.get("max_cpus_per_node", 0),
        min_cpus_per_node=node_cpus.get("min_cpus_per_node", 0),
        override=override,
    )


_SEARCH_CONCURRENCY_DEFAULT = 8  # W: polite default max concurrent array tasks


def derive_search_concurrency(
    override: int | None = None,
    qos_max_jobs: int | None = None,
) -> int:
    """Max concurrent search-array tasks (W) -- the array ``%W`` throttle.

    W is the SPEED lever: total search time ~= total_cost / (W * C). It is a
    polite *intent* cap, never tuned from a prior run -- the config default (8)
    or a positive ``override`` (from ``MultinodeConfig.array_throttle``). When a
    QOS ``MaxJobs`` limit is supplied it caps W (values < 1 mean "no limit" and
    are ignored); SLURM enforces QOS at runtime regardless. The further cap
    ``W <= T`` -- can't run more concurrently than there are tasks -- is applied
    where the array header is emitted, once T is known.
    """
    if override is not None:
        if override < 1:
            raise ValueError(
                f"array_throttle override must be >= 1, got {override}"
            )
        base = override
    else:
        base = _SEARCH_CONCURRENCY_DEFAULT
    if qos_max_jobs is not None and qos_max_jobs >= 1:
        base = min(base, qos_max_jobs)
    return max(1, base)


def resolve_search_concurrency(
    partition: str,
    override: int | None = None,
) -> int:
    """Resolve W for *partition*: config/override, capped by the discovered QOS.

    Thin submit-time wiring over :func:`derive_search_concurrency` and
    :func:`convgeno.slurm.discovery.detect_qos_max_jobs`. When the QOS cannot be
    determined the probe returns ``None`` and W falls back to the config
    default / override untouched.
    """
    return derive_search_concurrency(
        override=override,
        qos_max_jobs=detect_qos_max_jobs(partition),
    )


_SEARCH_WAVES_DEFAULT = 4  # K: waves of W tasks -> T = K * W buckets


def derive_search_waves(override: int | None = None) -> int:
    """Number of waves (K) the T = K*W buckets fan out over.

    K is a granularity/balancing knob: a larger K makes more, smaller buckets,
    giving LPT finer control and better backfill at the cost of more array
    tasks. Config default 4, or a positive ``override`` (from
    ``MultinodeConfig``). It is never capped by discovery -- K only multiplies
    the task count, not the concurrency (that is W).
    """
    if override is not None:
        if override < 1:
            raise ValueError(f"waves (K) override must be >= 1, got {override}")
        return override
    return _SEARCH_WAVES_DEFAULT


def derive_search_task_count(
    num_commands: int,
    concurrency: int,
    waves: int,
    max_array_size: int | None = None,
) -> int:
    """Number of LPT buckets / search-array tasks (T).

    ``T = K*W`` (waves x concurrency), clamped into ``[W, min(n^2,
    MaxArraySize)]``:

    * **upper bound (HARD)** -- T must not exceed the number of commands
      (``num_commands`` = n^2; more buckets would make empty tasks) nor SLURM's
      ``MaxArraySize`` (a larger array is rejected);
    * **lower bound (SOFT)** -- aim for at least W tasks so every concurrency
      slot can be filled. Because ``K >= 1`` this holds automatically whenever
      the ceiling allows; when the ceiling is below W (tiny dataset or a very
      small ``MaxArraySize``) the hard ceiling wins and W is capped to T where
      the array header is emitted.

    Hence ``T = min(K*W, num_commands, MaxArraySize?)``, clamped ``>= 1``. A
    ``max_array_size`` of ``None`` or ``< 1`` means "no array-size cap known".
    """
    if num_commands < 1:
        raise ValueError(f"num_commands must be >= 1, got {num_commands}")
    if concurrency < 1:
        raise ValueError(f"concurrency (W) must be >= 1, got {concurrency}")
    if waves < 1:
        raise ValueError(f"waves (K) must be >= 1, got {waves}")
    ceiling = num_commands
    if max_array_size is not None and max_array_size >= 1:
        ceiling = min(ceiling, max_array_size)
    return max(1, min(waves * concurrency, ceiling))


def resolve_search_task_count(
    num_commands: int,
    concurrency: int,
    waves: int,
) -> int:
    """Resolve T with ``MaxArraySize`` pulled from cluster discovery.

    Thin submit-time wiring over :func:`derive_search_task_count` and
    :func:`convgeno.slurm.discovery.detect_max_array_size`. When the array-size
    limit cannot be determined the probe returns ``None`` and only the
    ``num_commands`` ceiling applies.
    """
    return derive_search_task_count(
        num_commands=num_commands,
        concurrency=concurrency,
        waves=waves,
        max_array_size=detect_max_array_size(),
    )


# ------------------------------------------------------------
# Per-task --time / --mem sizing (from the biggest LPT bucket)
# ------------------------------------------------------------
# TIME. The bucket cost is in bytes^2 (sum of |S_i|*|S_j|). A task keeps all C
# cores busy -- C/p concurrent p-threaded diamonds -- so it burns cost at
# (C * throughput_const) bytes^2/second:
#
#     seconds = max_bucket_cost / (C * throughput_const) * margin
#
# throughput_const is a per-CORE diamond alignment rate (bytes^2 per
# core-second); p cancels (it repackages the same core-seconds, so it drives
# memory, not time). The default is a deliberately CONSERVATIVE placeholder:
# under-estimating time risks a task timeout (recoverable -- the resume
# completeness gate reports missing pairs + a resubmit line), while
# over-estimating only costs queue latency. CALIBRATE after the first real run:
#     throughput_const ~= max_bucket_cost / (C * observed_task_seconds) * margin
#
# MEMORY. C/p concurrent diamonds each load their target database, so peak RAM
# is concurrency * per_command + base. per_command uses the largest target DB
# (bytes) but is floored: for small proteomes (~10 MB) diamond's fixed working
# set (seed index + query block) dwarfs the loaded DB, so DB size alone would
# under-provision; db_safety only takes over for unusually large databases.

_SEARCH_THROUGHPUT_BYTES2_PER_CORE_SEC = 1.0e12  # conservative; CALIBRATE
_SEARCH_TIME_MARGIN = 1.5
_SEARCH_TIME_MIN_SEC = 900  # 15 min floor (diamond startup + DB load + I/O)
_SEARCH_TIME_MAX_SEC = 259200  # 72 h safety cap; hitting it => raise T

_SEARCH_MEM_DB_SAFETY = 4.0  # in-memory DB inflation over FASTA bytes
_SEARCH_MEM_PER_COMMAND_FLOOR_MB = 2048  # diamond per-process working set
_SEARCH_MEM_BASE_MB = 2048  # OS + orchestration overhead
_SEARCH_MEM_ROUND_MB = 1024
_SEARCH_MEM_FLOOR_MB = 4096


def within_task_concurrency(
    cpus_per_task: int,
    threads_per_command: int = 1,
) -> int:
    """Concurrent diamond commands per task = ``C // p`` (>= 1).

    Each emitted command uses ``p`` threads (``-p``, expected 1); running
    ``C // p`` of them at once keeps all C cores busy. ``p`` is detected from
    the emitted command later; the default 1 means concurrency == C.
    """
    if cpus_per_task < 1:
        raise ValueError(f"cpus_per_task must be >= 1, got {cpus_per_task}")
    if threads_per_command < 1:
        raise ValueError(
            f"threads_per_command must be >= 1, got {threads_per_command}"
        )
    return max(1, cpus_per_task // threads_per_command)


def derive_search_walltime(
    max_bucket_cost: int,
    cpus_per_task: int,
    throughput_const: float = _SEARCH_THROUGHPUT_BYTES2_PER_CORE_SEC,
    margin: float = _SEARCH_TIME_MARGIN,
) -> str | None:
    """Per-task walltime ``HH:MM:SS`` from the biggest bucket's cost.

    ``seconds = max_bucket_cost / (cpus_per_task * throughput_const) * margin``,
    clamped to ``[15 min, 72 h]``. ``throughput_const`` is bytes^2 per
    core-second (see the section comment; CALIBRATE after the first run).
    Returns ``None`` when not estimable (no cost / cores / rate) so the caller
    can fall back to the user default walltime.
    """
    if max_bucket_cost <= 0 or cpus_per_task <= 0 or throughput_const <= 0:
        return None
    seconds = max_bucket_cost / (cpus_per_task * throughput_const) * margin
    seconds = int(math.ceil(seconds))
    seconds = max(_SEARCH_TIME_MIN_SEC, min(_SEARCH_TIME_MAX_SEC, seconds))
    return _seconds_to_hms(seconds)


def derive_search_memory_mb(
    mem_determinant_bytes: int,
    concurrency: int,
    db_safety: float = _SEARCH_MEM_DB_SAFETY,
    per_command_floor_mb: int = _SEARCH_MEM_PER_COMMAND_FLOOR_MB,
    base_mb: int = _SEARCH_MEM_BASE_MB,
) -> int | None:
    """Per-task memory (MB) from the largest target database.

    ``mem = concurrency * max(largest_DB * db_safety, per_command_floor) +
    base``, rounded up to the nearest 1 GB with a 4 GB floor. ``concurrency``
    is ``C // p`` (concurrent diamonds, each loading its target DB). The
    per-command floor covers diamond's fixed working set, which dominates the
    small proteome DBs in this regime; ``db_safety`` only bites for unusually
    large databases. Returns ``None`` when not estimable.
    """
    if mem_determinant_bytes <= 0 or concurrency <= 0:
        return None
    largest_db_mb = mem_determinant_bytes / (1024 * 1024)
    per_command_mb = max(largest_db_mb * db_safety, per_command_floor_mb)
    raw = concurrency * per_command_mb + base_mb
    rounded = math.ceil(raw / _SEARCH_MEM_ROUND_MB) * _SEARCH_MEM_ROUND_MB
    return max(_SEARCH_MEM_FLOOR_MB, int(rounded))


# ------------------------------------------------------------
# Generate-time search-array sizing (one bundle: C, W, K, T, time, mem)
# ------------------------------------------------------------
# The search script's SBATCH header is static (written at generate time,
# before prepare runs), so its sizing is derived NOW from the INPUT FASTA
# sizes + cluster discovery -- exactly the inputs the plan calls for. T is
# deterministic (K*W clamped by n^2/MaxArraySize), so the same T is baked into
# the prepare script (for the manifest) and here (for --array): the manifest
# and the array width always agree. The --time/--mem values are ESTIMATES from
# input sizes (real Species*.fa sizes drive the actual manifest at prepare
# time, but the header only needs a sized-with-margin walltime/memory).

_SEARCH_FALLBACK_MEM = "16000M"  # off-cluster / no-size fallback only


def _input_fasta_sizes(input_dir: str) -> list[int]:
    """Byte sizes of the proteome FASTAs in ``input_dir`` (empty if unreadable)."""
    try:
        directory = Path(input_dir)
        if not directory.is_dir():
            return []
        return [
            entry.stat().st_size
            for entry in directory.iterdir()
            if entry.is_file() and entry.suffix.lower() in _FASTA_EXTENSIONS
        ]
    except OSError:
        return []


def build_all_pairs_search_commands(sizes: list[int]) -> list[SearchCommand]:
    """All ``n^2`` ordered species pairs as ``SearchCommand``s from proteome sizes.

    ``line_index`` is a synthetic enumeration (not tied to a real commands
    file); this feeds the generate-time LPT cost/mem ESTIMATE only, when the
    emitted commands don't exist yet. The real manifest is built at prepare
    time from ``_search_commands.txt`` via :func:`build_search_commands`.
    """
    commands: list[SearchCommand] = []
    idx = 0
    for query in sizes:
        for db in sizes:
            commands.append(
                SearchCommand(line_index=idx, query_bytes=query, db_bytes=db)
            )
            idx += 1
    return commands


@dataclass(frozen=True)
class SearchArraySizing:
    """Derived search-array parameters for one generated script."""

    cpus: int  # C -> --cpus-per-task
    concurrency: int  # W -> --array ...%W throttle
    waves: int  # K
    tasks: int  # T -> --array=0-(T-1)
    within_task_parallel: int  # C/p concurrent diamonds per task
    threads_per_command: int  # p
    time_limit: str  # per-task --time
    mem: str  # per-task --mem


def compute_search_array_sizing(
    input_sizes: list[int],
    max_cpus_per_node: int,
    min_cpus_per_node: int,
    qos_max_jobs: int | None,
    max_array_size: int | None,
    *,
    fallback_time: str,
    fallback_mem: str,
    overrides: MultinodeConfig | None = None,
) -> SearchArraySizing:
    """Pure derivation of the whole C/W/K/T + time/mem bundle.

    ``input_sizes`` are proteome FASTA byte sizes (empty off-cluster). All
    cluster numbers are passed in so this stays subprocess-free and testable;
    :func:`derive_search_array_sizing` is the discovery-backed wrapper.
    ``overrides`` (a :class:`MultinodeConfig`) pins any value the user set.
    """
    mn = overrides or MultinodeConfig()
    cpus = derive_search_cpus(
        max_cpus_per_node, min_cpus_per_node, override=mn.search_cpus
    )
    concurrency = derive_search_concurrency(
        override=mn.array_throttle, qos_max_jobs=qos_max_jobs
    )
    waves = derive_search_waves(override=mn.waves)
    threads_per_command = mn.threads_per_command or 1
    within = within_task_concurrency(cpus, threads_per_command)

    n = len(input_sizes)
    if n > 0:
        tasks = derive_search_task_count(
            n * n, concurrency, waves, max_array_size
        )
        lpt = lpt_partition(build_all_pairs_search_commands(input_sizes), tasks)
        max_bucket_cost = lpt.max_bucket_cost
        mem_determinant = lpt.mem_determinant_bytes
    else:
        # Off-cluster / no data: still size the array, fall back on time/mem.
        base = waves * concurrency
        tasks = min(base, max_array_size) if max_array_size else base
        tasks = max(1, tasks)
        max_bucket_cost = 0
        mem_determinant = 0

    throughput = mn.throughput_const or _SEARCH_THROUGHPUT_BYTES2_PER_CORE_SEC
    margin = mn.time_margin or _SEARCH_TIME_MARGIN

    if mn.search_time_limit:
        time_limit = mn.search_time_limit
    else:
        time_limit = (
            derive_search_walltime(max_bucket_cost, cpus, throughput, margin)
            or fallback_time
        )

    if mn.search_mem:
        mem = mn.search_mem
    else:
        mem_mb = derive_search_memory_mb(mem_determinant, within)
        mem = f"{mem_mb}M" if mem_mb is not None else fallback_mem

    return SearchArraySizing(
        cpus=cpus,
        concurrency=concurrency,
        waves=waves,
        tasks=tasks,
        within_task_parallel=within,
        threads_per_command=threads_per_command,
        time_limit=time_limit,
        mem=mem,
    )


def derive_search_array_sizing(config: PipelineConfig) -> SearchArraySizing:
    """Discovery-backed :func:`compute_search_array_sizing` for *config*.

    Probes the partition (node CPUs, QOS ``MaxJobs``, ``MaxArraySize``) and
    reads the input FASTA sizes, then defers to the pure computation. All
    probes degrade to safe fallbacks off-cluster.
    """
    node = detect_node_cpus(config.slurm.partition)
    input_sizes = (
        _input_fasta_sizes(config.orthofinder.input_dir)
        if config.orthofinder is not None
        else []
    )
    return compute_search_array_sizing(
        input_sizes,
        max_cpus_per_node=node.get("max_cpus_per_node", 0),
        min_cpus_per_node=node.get("min_cpus_per_node", 0),
        qos_max_jobs=detect_qos_max_jobs(config.slurm.partition),
        max_array_size=detect_max_array_size(),
        fallback_time=config.slurm.time_limit,
        fallback_mem=config.slurm.mem or _SEARCH_FALLBACK_MEM,
        overrides=config.multinode,
    )


# Reusable bash helpers embedded in the search-array script. Kept as raw
# module constants (not f-strings) so their shell ``{}``/``$`` need no escaping.

# Resolve a writable scratch base (self-heal owner-write, then write-probe).
# Lifted from the single-node generator so both paths behave identically.
_RESOLVE_SCRATCH_FUNC = r"""
resolve_scratch_base() {
    local base="$1"
    if [ -d "$base" ] && [ ! -w "$base" ] && [ -O "$base" ]; then
        chmod u+w "$base" 2>/dev/null || true
    fi
    mkdir -p "$base" 2>/dev/null || true
    if [ -d "$base" ]; then
        local probe="$base/.convgeno_wtest.$$"
        if ( : > "$probe" ) 2>/dev/null; then
            rm -f "$probe"
            printf '%s' "$base"
            return 0
        fi
    fi
    return 1
}
"""

# Per-command worker (run concurrently via xargs -P). Idempotent skip on a
# valid existing Blast{i}_{j}.txt.gz in the shared WorkingDirectory; otherwise
# run the emitted diamond command with its -o directory redirected to
# $OUTPUT_BASE (scratch or the WorkingDirectory) -- inputs (-q/-d) stay shared.
# Never returns non-zero (so xargs runs the whole bucket); outcome is recorded
# as a marker file the parent counts.
_RUN_ONE_FUNC = r"""
run_one() {
    local CMD="$1"
    local PAIR
    if [[ "$CMD" =~ Blast([0-9]+)_([0-9]+) ]]; then
        PAIR="${BASH_REMATCH[0]}"
    else
        : > "$STATUS_DIR/failed/noblast_$$_$RANDOM"
        return 0
    fi
    local FINAL="$WORK_DIR/${PAIR}.txt.gz"
    if [ -s "$FINAL" ] && gzip -t "$FINAL" 2>/dev/null; then
        : > "$STATUS_DIR/skipped/$PAIR"
        return 0
    fi
    local RUN
    RUN="$(printf '%s' "$CMD" | sed -E "s#(-o[[:space:]]+)[^[:space:]]*/([^/[:space:]]+)#\1$OUTPUT_BASE/\2#")"
    if eval "$RUN" >/dev/null 2>&1; then
        : > "$STATUS_DIR/ran/$PAIR"
    else
        : > "$STATUS_DIR/failed/$PAIR"
    fi
    return 0
}
export -f run_one
"""


def generate_search_array_script(
    config: PipelineConfig,
    runtime: CondaRuntimeConfig | None = None,
    sizing: SearchArraySizing | None = None,
) -> str:
    """Generate the SLURM array job for the OrthoFinder search phase (v2).

    A backfill-friendly array of ``T`` LPT-balanced tasks (``--array=0-(T-1)%W``,
    ``--cpus-per-task=C``). Each task reads its manifest
    (``search_task_<id>.txt``) and runs ``C/p`` raw ``diamond blastp`` commands
    concurrently, idempotently skipping any valid existing ``Blast{i}_{j}.txt.gz``.

    Outputs are written to scratch (when ``slurm.scratch_dir`` is set: preferred
    scratch -> ``/tmp`` -> **direct-to-shared** if neither is writable) and
    rsynced back to the WorkingDirectory, with a SIGTERM/SIGINT salvage trap and
    ephemeral/persistent cleanup -- the single-node scratch pattern. Sizing is
    derived once (or passed in for consistency with the prepare script).
    """
    _validate_has_orthofinder(config)
    assert config.orthofinder is not None
    resolved_runtime = _resolve_runtime(config, runtime)
    if sizing is None:
        sizing = derive_search_array_sizing(config)

    throttle = min(sizing.concurrency, sizing.tasks)
    array_spec = f"0-{sizing.tasks - 1}%{throttle}"
    array_output = config.slurm.output_pattern.replace("%j", "%A_%a")
    array_error = config.slurm.error_pattern.replace("%j", "%A_%a")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header = _build_sbatch_header(
        config,
        overrides={
            "--job-name": "convgeno_of_search",
            "--array": array_spec,
            "--nodes": "1",
            "--ntasks": "1",
            "--cpus-per-task": str(sizing.cpus),
            "--time": sizing.time_limit,
            "--mem": sizing.mem,
            "--mem-per-cpu": None,
            "--output": array_output,
            "--error": array_error,
        },
    )
    bootstrap_block = render_conda_bootstrap(resolved_runtime)

    scratch_dir = config.slurm.scratch_dir
    if scratch_dir is not None:
        scratch_setup_block = f"""
# ---- Scratch: write Blast outputs locally, rsync back to shared ----
PREFERRED_SCRATCH="{scratch_dir}"
FALLBACK_SCRATCH="/tmp/scratch"
{_RESOLVE_SCRATCH_FUNC}
SCRATCH_BASE="$(resolve_scratch_base "$PREFERRED_SCRATCH")"
if [ -z "$SCRATCH_BASE" ]; then
    echo "WARNING: preferred scratch '$PREFERRED_SCRATCH' not writable; trying '$FALLBACK_SCRATCH'"
    SCRATCH_BASE="$(resolve_scratch_base "$FALLBACK_SCRATCH")"
fi
if [ -n "$SCRATCH_BASE" ]; then
    USE_SCRATCH=1
    JOB_SCRATCH="$SCRATCH_BASE/${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"
    OUTPUT_BASE="$JOB_SCRATCH/blast_out"
    export TMPDIR="$JOB_SCRATCH/tmp"
    mkdir -p "$OUTPUT_BASE" "$TMPDIR"
    echo "Scratch base: $SCRATCH_BASE (Blast outputs -> $OUTPUT_BASE)"
    salvage_scratch() {{
        trap - SIGTERM SIGINT
        echo "Interrupted — salvaging Blast outputs from scratch to WorkingDirectory..."
        rsync -a "$OUTPUT_BASE"/ "$WORK_DIR"/ 2>/dev/null || true
    }}
    trap salvage_scratch SIGTERM SIGINT
else
    USE_SCRATCH=0
    OUTPUT_BASE="$WORK_DIR"
    echo "No writable scratch; writing Blast outputs directly to the shared WorkingDirectory."
fi
"""
        if config.slurm.is_ephemeral_scratch:
            scratch_cleanup_block = """
if [ "$USE_SCRATCH" = "1" ] && [ -d "$JOB_SCRATCH" ]; then
    echo "Cleaning up ephemeral scratch: $JOB_SCRATCH"
    rm -rf "$JOB_SCRATCH" || echo "WARNING: could not remove $JOB_SCRATCH (node reclaims it)"
fi
"""
        else:
            scratch_cleanup_block = """
if [ "$USE_SCRATCH" = "1" ]; then
    echo "Persistent scratch left in $JOB_SCRATCH (cluster purge policy reclaims it)."
fi
"""
    else:
        scratch_setup_block = """
# ---- Scratch not configured: write directly to the shared WorkingDirectory ----
USE_SCRATCH=0
OUTPUT_BASE="$WORK_DIR"
echo "Scratch not configured; writing Blast outputs directly to the shared WorkingDirectory."
"""
        scratch_cleanup_block = ""

    script = f"""\
#!/bin/bash
# ============================================================
# OrthoFinder search array (raw diamond blastp) — multi-node mode
# Generated by convgeno on {timestamp}
# ============================================================

{header}

set -uo pipefail

{bootstrap_block}

echo "============================================================"
echo "convgeno: OrthoFinder search array task started"
echo "Array job: ${{SLURM_ARRAY_JOB_ID:-?}}  Task: ${{SLURM_ARRAY_TASK_ID:-?}}"
echo "Node:      $HOSTNAME"
echo "Start:     $(date)"
echo "============================================================"

START_SECONDS=$SECONDS

OUTPUT_DIR="{config.orthofinder.output_dir}"
OUTPUT_PARENT="$(dirname "$OUTPUT_DIR")"
RUN_NAME="$(basename "$OUTPUT_DIR")"
WORK_DIR_FILE="$OUTPUT_PARENT/${{RUN_NAME}}_working_dir_path.txt"
MANIFEST_DIR="$OUTPUT_PARENT/${{RUN_NAME}}{_SEARCH_MANIFEST_DIRNAME_SUFFIX}"

if [ ! -f "$WORK_DIR_FILE" ]; then
    echo "ERROR: WorkingDirectory pointer not found: $WORK_DIR_FILE (prepare failed?)" >&2
    exit 1
fi
WORK_DIR="$(cat "$WORK_DIR_FILE")"
if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: WorkingDirectory not found or unreadable: $WORK_DIR" >&2
    exit 1
fi

TASK_MANIFEST="$MANIFEST_DIR/search_task_${{SLURM_ARRAY_TASK_ID}}.txt"
if [ ! -f "$TASK_MANIFEST" ]; then
    echo "ERROR: task manifest not found: $TASK_MANIFEST" >&2
    echo "The prepare job builds these; check that it completed." >&2
    exit 1
fi
{scratch_setup_block}
STATUS_DIR="$(mktemp -d)"
mkdir -p "$STATUS_DIR/ran" "$STATUS_DIR/skipped" "$STATUS_DIR/failed"

{_RUN_ONE_FUNC}
export WORK_DIR OUTPUT_BASE STATUS_DIR

CONCURRENCY={sizing.within_task_parallel}
N_CMDS=$(grep -c . "$TASK_MANIFEST" || true)
echo "Task ${{SLURM_ARRAY_TASK_ID}}: $N_CMDS commands, $CONCURRENCY concurrent (C/p), writing to $OUTPUT_BASE"

# Run the bucket: C/p diamonds at a time, one command per worker (NUL-delimited
# so command spaces are preserved). Per-command outcomes are marker files.
tr '\\n' '\\0' < "$TASK_MANIFEST" | xargs -0 -r -n1 -P "$CONCURRENCY" bash -c 'run_one "$1"' _

if [ "$USE_SCRATCH" = "1" ]; then
    echo "Copying Blast outputs: scratch -> WorkingDirectory"
    rsync -a "$OUTPUT_BASE"/ "$WORK_DIR"/
fi

RAN=$(find "$STATUS_DIR/ran" -type f | wc -l)
SKIPPED=$(find "$STATUS_DIR/skipped" -type f | wc -l)
FAILED=$(find "$STATUS_DIR/failed" -type f | wc -l)

ELAPSED=$(( SECONDS - START_SECONDS ))
echo "============================================================"
echo "convgeno: search task ${{SLURM_ARRAY_TASK_ID}} finished"
echo "ran=$RAN skipped=$SKIPPED failed=$FAILED  (of $N_CMDS commands)"
echo "Duration: ${{ELAPSED}}s   End: $(date)"
echo "============================================================"

rm -rf "$STATUS_DIR"
{scratch_cleanup_block}
if [ "$FAILED" -gt 0 ]; then
    echo "ERROR: $FAILED search command(s) failed in task ${{SLURM_ARRAY_TASK_ID}}." >&2
    echo "Fix the cause and resubmit this array id; completed Blast files are skipped." >&2
    exit 1
fi
exit 0
"""
    return script


def generate_resume_script(
    config: PipelineConfig,
    runtime: CondaRuntimeConfig | None = None,
) -> str:
    """Generate SLURM script for OrthoFinder's resume phase (``-b``).

    Runs on a single node after all search array tasks complete. Resumes
    OrthoFinder from pre-computed search results, performing clustering,
    MSA, and tree inference.

    The generated script:

    * Sets ``TMPDIR`` to a project-local temp directory to avoid filling
      ``/dev/shm`` or ``/tmp``.
    * Cleans up ``$TMPDIR`` on success, preserves it on failure.
    * Uses ``orthofinder.analysis_threads`` for ``-a`` (falling back to 1
      when unset).

    NOTE: open-file-limit handling was removed from this script; the fd
    feasibility gate is reintroduced as a separate step.
    """
    _validate_has_orthofinder(config)
    assert config.orthofinder is not None
    resolved_runtime = _resolve_runtime(config, runtime)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header = _build_sbatch_header(
        config,
        overrides={
            "--job-name": "convgeno_of_resume",
            "--nodes": "1",
            "--ntasks": "1",
            "--cpus-per-task": str(config.slurm.cpus_per_task),
            "--time": config.slurm.time_limit,
        },
    )
    bootstrap_block = render_conda_bootstrap(resolved_runtime)

    # Resolve analysis threads (-a): use the configured value, falling back
    # to 1 when unset. The open-file-limit-driven auto-heuristic was removed;
    # the fd feasibility gate (a later step) will own fd sizing.
    if config.orthofinder.analysis_threads is not None:
        analysis_threads = config.orthofinder.analysis_threads
    else:
        analysis_threads = 1

    extra = " ".join(config.orthofinder.extra_args)
    extra_suffix = f" {extra}" if extra else ""

    script = f"""\
#!/bin/bash
# ============================================================
# OrthoFinder resume phase (-b) — multi-node mode
# Generated by convgeno on {timestamp}
# ============================================================

{header}

set -euo pipefail

{bootstrap_block}

echo "============================================================"
echo "convgeno: OrthoFinder resume phase started"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $HOSTNAME"
echo "Start:     $(date)"
echo "============================================================"

START_SECONDS=$SECONDS

OUTPUT_DIR="{config.orthofinder.output_dir}"
OUTPUT_PARENT="$(dirname "$OUTPUT_DIR")"
RUN_NAME="$(basename "$OUTPUT_DIR")"
WORK_DIR_FILE="$OUTPUT_PARENT/${{RUN_NAME}}_working_dir_path.txt"

# The prepare script wrote the WorkingDirectory path to a pointer file
# in $OUTPUT_PARENT. Read it directly — no find-fallback, since a
# missing pointer means the prepare job failed and we should abort.
if [ ! -f "$WORK_DIR_FILE" ]; then
    echo "ERROR: WorkingDirectory pointer not found: $WORK_DIR_FILE"
    echo "The prepare job may have failed."
    exit 1
fi

WORK_DIR="$(cat "$WORK_DIR_FILE")"

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: WorkingDirectory not found or unreadable: $WORK_DIR"
    exit 1
fi

# ---- Project-local TMPDIR ----
export TMPDIR="$OUTPUT_PARENT/${{RUN_NAME}}_tmp"
mkdir -p "$TMPDIR"
echo "TMPDIR=$TMPDIR"

TOTAL_THREADS={config.orthofinder.search_threads}
ANALYSIS_THREADS={analysis_threads}

echo "Resuming OrthoFinder from: $WORK_DIR"
echo "Search threads (-t): $TOTAL_THREADS"
echo "Analysis threads (-a): $ANALYSIS_THREADS"
echo ""

orthofinder -b "$WORK_DIR" -t "$TOTAL_THREADS" -a "$ANALYSIS_THREADS"{extra_suffix}

OF_EXIT=$?

ELAPSED=$(( SECONDS - START_SECONDS ))
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "============================================================"
echo "convgeno: resume phase finished"
echo "Exit code: $OF_EXIT"
echo "End:       $(date)"
echo "Duration:  ${{HOURS}}h ${{MINUTES}}m ${{SECS}}s"
echo "============================================================"

if [ "$OF_EXIT" -eq 0 ]; then
    echo "OrthoFinder completed successfully. Cleaning up TMPDIR: $TMPDIR"
    rm -rf "$TMPDIR"

    if [ ! -f "$OUTPUT_DIR"/*/Orthogroups/Orthogroups.tsv ] 2>/dev/null; then
        echo "WARNING: Expected output file Orthogroups.tsv not found in $OUTPUT_DIR"
        echo "Check $OUTPUT_DIR for a Results_* directory."
    fi

    echo "OrthoFinder multi-node run completed successfully."
else
    echo "OrthoFinder exited with code $OF_EXIT. Preserving TMPDIR for debugging: $TMPDIR"
fi

exit "$OF_EXIT"
"""
    return script
