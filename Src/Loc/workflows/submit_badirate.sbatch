#!/usr/bin/env bash
#SBATCH --job-name=badirate
#SBATCH --output=Data/logs/badirate/badirate_%j.out
#SBATCH --error=Data/logs/badirate/badirate_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -euo pipefail

CONFIG="Src/Loc/configs/example_config.yaml"
TOOLS="Src/Loc/configs/tool_paths.yaml"

python Src/Loc/scripts/run_badirate.py --config "$CONFIG" --tools "$TOOLS"
