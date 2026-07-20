#!/usr/bin/env bash
#SBATCH --job-name=orthofinder
#SBATCH --output=Data/logs/orthofinder/orthofinder_%j.out
#SBATCH --error=Data/logs/orthofinder/orthofinder_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#
# submit_orthofinder.sbatch — Run OrthoFinder on SLURM.
# Adjust --time, --mem, --cpus-per-task based on species count.

set -euo pipefail

CONFIG="Src/Loc/configs/example_config.yaml"
TOOLS="Src/Loc/configs/tool_paths.yaml"

python Src/Loc/scripts/run_orthofinder.py --config "$CONFIG" --tools "$TOOLS"

# Prepare inputs for downstream tools
python Src/Loc/scripts/prepare_cafe_inputs.py --config "$CONFIG"
