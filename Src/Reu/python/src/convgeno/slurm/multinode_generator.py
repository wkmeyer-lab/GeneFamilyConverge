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

from convgeno.slurm.config import PipelineConfig, normalize_optional_account
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
) -> str:
    """Generate SLURM script for OrthoFinder's prepare phase (``-op``).

    Runs on a single node. Captures the DIAMOND/BLAST commands to a file
    for the array job.
    """
    _validate_has_orthofinder(config)
    assert config.orthofinder is not None  # for type checker
    resolved_runtime = _resolve_runtime(config, runtime)

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


def generate_search_array_script(
    config: PipelineConfig,
    runtime: CondaRuntimeConfig | None = None,
) -> str:
    """Generate the SLURM array job for the OrthoFinder search phase.

    The v1 search array has been removed. The v2 implementation (LPT-balanced
    buckets, derived C/W/K/T sizing, within-task ``C/p`` concurrency, idempotent
    ``Blast`` skip) is not yet built — see the plan's "Search-array sizing"
    section.
    """
    raise NotImplementedError(
        "Multi-node search array is being reimplemented (v2: LPT buckets + "
        "derived C/W/K/T sizing). It is not yet available."
    )


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
