#!/usr/bin/env bash
# run_local.sh — Execute the full pipeline locally (Mac/Linux workstation).
#
# Usage:
#   bash Src/Loc/workflows/run_local.sh Src/Loc/configs/example_config.yaml
#
# This is the manual execution path. Once an orchestrator (Snakemake/Nextflow)
# is chosen, this script becomes a reference for the step ordering.

set -euo pipefail

CONFIG="${1:?Usage: run_local.sh <config.yaml>}"
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

echo "=== Running OrthoFinder ==="
python Src/Loc/scripts/run_orthofinder.py --config "$CONFIG" --tools "$TOOLS"

echo "=== Making species tree ultrametric (r8s) ==="
python Src/Loc/scripts/make_tree_ultrametric.py --config "$CONFIG" --tools "$TOOLS"

echo "=== Preparing CAFE-5 inputs ==="
python Src/Loc/scripts/prepare_cafe_inputs.py --config "$CONFIG"

echo "=== Building categorical phenotype tree (CAFE -y) ==="
Rscript Src/Loc/scripts/make_categorical_phenotype_tree.R --config "$CONFIG"

echo "=== Running CAFE-5 ==="
python Src/Loc/scripts/run_cafe.py --config "$CONFIG" --tools "$TOOLS"

echo "=== Summarizing results (R) ==="
Rscript Src/Loc/scripts/summarize_results.R --config "$CONFIG"

echo "=== Generating figures (R) ==="
Rscript Src/Loc/scripts/make_figures.R --config "$CONFIG"

echo "=== Pipeline complete ==="
