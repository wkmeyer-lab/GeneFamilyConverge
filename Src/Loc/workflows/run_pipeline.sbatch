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
# the heavy tools (OrthoFinder, CAFE-5) as separate jobs
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

echo "=== Filtering proteomes to longest isoform per gene ==="
# filter_isoforms.py takes explicit paths, not --config: keep these in sync with
# inputs.proteome_dir and orthofinder.input_dir in $CONFIG; --format matches
# inputs.source_format.
python Src/Loc/scripts/filter_isoforms.py dir \
    Data/raw/proteomes \
    Data/interim/cleaned_proteomes \
    --format ensembl \
    --stats-json Data/logs/isoform_filter_stats.json

echo "=== Submitting OrthoFinder ==="
# The OrthoFinder job also makes the species tree ultrametric (r8s) and
# prepares the CAFE-5 inputs once OrthoFinder finishes (see submit_orthofinder.sbatch).
JOB_OF=$(sbatch --parsable Src/Loc/workflows/submit_orthofinder.sbatch)
echo "OrthoFinder job: $JOB_OF"

echo "=== Submitting CAFE-5 (after OrthoFinder) ==="
JOB_CAFE=$(sbatch --parsable --dependency=afterok:${JOB_OF} Src/Loc/workflows/submit_cafe.sbatch)
echo "CAFE-5 job: $JOB_CAFE"

echo "=== Pipeline submitted. CAFE-5 runs after OrthoFinder. ==="
echo "After it completes, run the R summarization step manually:"
echo "  Rscript Src/Loc/scripts/summarize_results.R --config $CONFIG"
