#!/usr/bin/env bash
#SBATCH --job-name=genefam-pipeline
#SBATCH --output=Data/logs/pipeline/pipeline_%j.out
#SBATCH --error=Data/logs/pipeline/pipeline_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#
# run_pipeline.sbatch — Orchestrate the full pipeline on Lehigh SLURM.
#
# This lightweight job runs the fast Python steps directly and submits
# the heavy tools (OrthoFinder, CAFE-5, BadiRate) as separate jobs
# with --dependency so they run in order.
#
# Usage:
#   sbatch Src/Loc/workflows/run_pipeline.sbatch
#
# TODO: Replace with Snakemake/Nextflow once orchestrator is chosen.
#       This script exists so you can run on HPC before that decision.

set -euo pipefail

CONFIG="Src/Loc/configs/example_config.yaml"
TOOLS="Src/Loc/configs/tool_paths.yaml"

echo "=== Validating inputs ==="
python Src/Loc/scripts/validate_inputs.py --config "$CONFIG" --tools "$TOOLS"

echo "=== Preparing proteomes ==="
python Src/Loc/scripts/prepare_proteomes.py --config "$CONFIG"

echo "=== Submitting OrthoFinder ==="
JOB_OF=$(sbatch --parsable Src/Loc/workflows/submit_orthofinder.sbatch)
echo "OrthoFinder job: $JOB_OF"

echo "=== Submitting CAFE-5 (after OrthoFinder) ==="
JOB_CAFE=$(sbatch --parsable --dependency=afterok:${JOB_OF} Src/Loc/workflows/submit_cafe.sbatch)
echo "CAFE-5 job: $JOB_CAFE"

echo "=== Submitting BadiRate (after OrthoFinder) ==="
JOB_BR=$(sbatch --parsable --dependency=afterok:${JOB_OF} Src/Loc/workflows/submit_badirate.sbatch)
echo "BadiRate job: $JOB_BR"

echo "=== Pipeline submitted. CAFE-5 and BadiRate run in parallel after OrthoFinder. ==="
echo "After both complete, run the R summarization step manually:"
echo "  Rscript Src/Loc/scripts/summarize_results.R --config $CONFIG"
