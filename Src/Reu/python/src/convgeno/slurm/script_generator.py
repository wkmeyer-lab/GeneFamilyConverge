"""SLURM batch script generation for external bioinformatics tools."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from convgeno.slurm.config import PipelineConfig
from convgeno.slurm.runtime import CondaRuntimeConfig, render_conda_bootstrap


def _orthofinder_command_with_paths(
    args: list[str],
    input_arg: str,
    output_arg: str,
) -> str:
    rebuilt: list[str] = []
    replace_next: str | None = None
    for arg in args:
        if replace_next == "-f":
            rebuilt.append(input_arg)
            replace_next = None
            continue
        if replace_next == "-o":
            rebuilt.append(output_arg)
            replace_next = None
            continue
        rebuilt.append(arg)
        if arg in {"-f", "-o"}:
            replace_next = arg
    return "orthofinder " + " ".join(rebuilt)


def generate_orthofinder_script(
    config: PipelineConfig,
    runtime: CondaRuntimeConfig | None = None,
) -> str:
    """Generate a SLURM batch script to run OrthoFinder.

    Parameters
    ----------
    config
        Pipeline configuration containing SLURM resources and
        OrthoFinder settings.
    runtime
        Conda runtime configuration for compute-node bootstrap.
        If ``None``, falls back to ``config.runtime``. At least one
        must be provided.

    Returns
    -------
    str
        Complete script ready to be written to disk and submitted
        with ``sbatch``.

    Raises
    ------
    ValueError
        If ``config.orthofinder`` is ``None`` or no runtime config
        is available.
    """
    if config.orthofinder is None:
        raise ValueError(
            "OrthoFinder settings not found in pipeline config. "
            "Add an 'orthofinder' section to pipeline_config.yaml."
        )

    resolved_runtime = runtime if runtime is not None else config.runtime
    if resolved_runtime is None:
        raise ValueError(
            "No runtime configuration provided. "
            "Run 'convgeno init' to detect and store conda paths."
        )

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sbatch_lines = config.slurm.to_sbatch_lines()
    sbatch_lines.insert(0, "#SBATCH --job-name=convgeno_orthofinder")
    sbatch_lines.append("#SBATCH --export=ALL")

    of_args = config.orthofinder.to_command_args()
    scratch_dir = config.slurm.scratch_dir
    scratch_enabled = scratch_dir is not None
    # ## NEW: Build scratch-aware OrthoFinder paths while preserving direct runs.
    if scratch_enabled:
        of_command = _orthofinder_command_with_paths(
            of_args,
            '"$EFFECTIVE_INPUT"',
            '"$EFFECTIVE_OUTPUT"',
        )
    else:
        of_command = _orthofinder_command_with_paths(
            of_args,
            config.orthofinder.input_dir,
            config.orthofinder.output_dir,
        )
    of_command_for_echo = of_command.replace('"', '\\"')

    scratch_setup_block = ""
    scratch_copy_inputs_block = ""
    scratch_rsync_back_block = ""
    scratch_cleanup_block = ""
    if scratch_enabled:
        scratch_setup_block = f"""
# ## NEW: Scratch space setup
JOB_SCRATCH="{scratch_dir}/${{SLURM_JOB_ID}}"
SCRATCH_INPUT="$JOB_SCRATCH/input_fastas"
SCRATCH_OUTPUT="$JOB_SCRATCH/orthofinder_output"
SCRATCH_TMP="$JOB_SCRATCH/tmp"
EFFECTIVE_INPUT="$SCRATCH_INPUT"
EFFECTIVE_OUTPUT="$SCRATCH_OUTPUT"

mkdir -p "$SCRATCH_INPUT" "$SCRATCH_TMP"

export TMPDIR="$SCRATCH_TMP"
export TEMP="$SCRATCH_TMP"
export TMP="$SCRATCH_TMP"

echo "Using scratch directory: $JOB_SCRATCH"
df -h "$JOB_SCRATCH" 2>/dev/null || true

# ## NEW: Trap to copy partial results on job cancellation or timeout.
cleanup_scratch() {{
    echo "Signal received — copying partial results from scratch..."
    if [ -d "$SCRATCH_OUTPUT" ]; then
        rsync -a --info=progress2 "$SCRATCH_OUTPUT"/ "$OUTPUT_DIR"/ 2>/dev/null || true
        echo "Partial results copied to: $OUTPUT_DIR"
    fi
}}
trap cleanup_scratch SIGTERM SIGINT EXIT
"""
        scratch_copy_inputs_block = """
# ## NEW: Copy validated FASTA inputs to scratch before running OrthoFinder.
echo "Copying FASTA files to scratch at: $(date)"
rsync -a --include="*.fa" --include="*.fasta" --include="*.faa" --exclude="*" "$INPUT_DIR"/ "$SCRATCH_INPUT"/
SCRATCH_FASTA_COUNT=$(find "$SCRATCH_INPUT" -maxdepth 1 -name "*.fa" -o -name "*.fasta" -o -name "*.faa" | wc -l)
echo "Copied $SCRATCH_FASTA_COUNT FASTA files to scratch"
echo "Scratch input size: $(du -sh "$SCRATCH_INPUT" | cut -f1)"
"""
        scratch_rsync_back_block = """
# ## NEW: Copy scratch results back before post-run validation.
if [ $EXIT_CODE -eq 0 ]; then
    echo "Copying results from scratch to final output at: $(date)"
    mkdir -p "$OUTPUT_DIR"
    rsync -a --info=progress2 "$SCRATCH_OUTPUT"/ "$OUTPUT_DIR"/
    echo "Results copied to: $OUTPUT_DIR"
fi
"""
        scratch_cleanup_block = """
# ## NEW: Clean up scratch after a successful result copy.
# Clean up scratch
if [ -d "$JOB_SCRATCH" ]; then
    echo "Cleaning up scratch directory: $JOB_SCRATCH"
    rm -rf "$JOB_SCRATCH"
fi
"""

    sbatch_block = "\n".join(sbatch_lines)
    bootstrap_block = render_conda_bootstrap(resolved_runtime)

    script = f"""\
#!/bin/bash
# ============================================================
# OrthoFinder SLURM batch script
# Generated by convgeno on {timestamp}
# ============================================================
#
# To submit:  sbatch this_script.sh
# To inspect: cat this_script.sh
#

{sbatch_block}

# ------------------------------------------------------------
# Environment setup
# ------------------------------------------------------------
set -euo pipefail

{bootstrap_block}

# Record start time for duration calculation
START_SECONDS=$SECONDS

# Verify OrthoFinder is available
if ! command -v orthofinder &> /dev/null; then
    echo "ERROR: orthofinder not found in conda environment"
    echo "Install with: conda install -c bioconda orthofinder"
    exit 1
fi

echo "OrthoFinder version:"
orthofinder -h 2>&1 | head -n 2 || true
echo ""

INPUT_DIR="{config.orthofinder.input_dir}"
OUTPUT_DIR="{config.orthofinder.output_dir}"
{scratch_setup_block}
# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

FASTA_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fa" -o -name "*.fasta" -o -name "*.faa" | wc -l)
if [ "$FASTA_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTA files (.fa, .fasta, .faa) found in $INPUT_DIR"
    exit 1
fi

echo "Input directory: $INPUT_DIR"
echo "FASTA files found: $FASTA_COUNT"
echo "Output directory: $OUTPUT_DIR"
echo ""

# OrthoFinder refuses to run when the non-default -o output directory
# already exists. Create only the PARENT directory, and fail loudly if
# the target output directory itself is present.
mkdir -p "$(dirname "$OUTPUT_DIR")"

if [ -e "$OUTPUT_DIR" ]; then
    echo "ERROR: OrthoFinder output directory already exists: $OUTPUT_DIR"
    echo "Choose a fresh output_dir in pipeline_config.yaml."
    exit 1
fi
{scratch_copy_inputs_block}

# ------------------------------------------------------------
# Run OrthoFinder
# ------------------------------------------------------------
echo "Running OrthoFinder..."
echo "Command: {of_command_for_echo}"
echo ""

{of_command}

EXIT_CODE=$?
{scratch_rsync_back_block}

# ------------------------------------------------------------
# Post-run validation
# ------------------------------------------------------------
ELAPSED=$(( SECONDS - START_SECONDS ))
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "============================================================"
echo "convgeno: OrthoFinder job finished"
echo "Exit code: $EXIT_CODE"
echo "End:       $(date)"
echo "Duration:  ${{HOURS}}h ${{MINUTES}}m ${{SECS}}s"
echo "============================================================"

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: OrthoFinder exited with code $EXIT_CODE"
    exit $EXIT_CODE
fi

# Check that key output files exist
if [ ! -f "$OUTPUT_DIR"/*/Orthogroups/Orthogroups.tsv ] 2>/dev/null; then
    echo "WARNING: Expected output file Orthogroups.tsv not found in $OUTPUT_DIR"
    echo "OrthoFinder may have written results to a subdirectory."
    echo "Check $OUTPUT_DIR for a Results_* directory."
fi

echo "OrthoFinder completed successfully."
{scratch_cleanup_block}
exit 0
"""
    return script


def write_script(content: str, path: Path | str) -> Path:
    """Write a generated SLURM script to disk and make it executable.

    Returns the resolved path.
    """
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | 0o755)
    return path
