#!/usr/bin/env bash
#SBATCH --job-name=cafe5
#SBATCH --output=Data/logs/cafe/cafe_%j.out
#SBATCH --error=Data/logs/cafe/cafe_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -euo pipefail

CONFIG="Src/Loc/configs/example_config.yaml"
TOOLS="Src/Loc/configs/tool_paths.yaml"

# Categorical phenotype tree (the CAFE -y multi-lambda tree). Lightweight R/ASR
# step (~1 core, seconds); its input cafe_input/species_tree.nwk is produced by
# the OrthoFinder job this one depends on (afterok). Requires R on PATH (e.g.
# `module load R`) with the tools/r-deps packages (ape, castor, expm, yaml).
Rscript Src/Loc/scripts/make_categorical_phenotype_tree.R --config "$CONFIG"

python Src/Loc/scripts/run_cafe.py --config "$CONFIG" --tools "$TOOLS"
