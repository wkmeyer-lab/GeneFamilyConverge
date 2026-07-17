#!/usr/bin/env bash
# SLURM wrapper for tools/throughput-constant/throughput_calibration.sh.
#
# This is a TEMPLATE. Edit the lines marked "# EDIT:" for your cluster, then:
#   sbatch tools/throughput-constant/throughput_calibration.sbatch
# Submit from the repository root (so SLURM_SUBMIT_DIR is the project root and
# pipeline_config.yaml is found).
#
# It reads the prepare job's outputs from pipeline_config.yaml, so run the
# prepare phase first:  convgeno orthofinder generate  &&  sbatch \
#   slurm_scripts/orthofinder_prepare.sh   (wait for COMPLETED).

# EDIT: match your partition and a single-core allocation. One core is the
# point -- we measure isolated single-thread throughput. Give ~1 node's
# per-core memory and enough time (100 sequential searches can take many
# hours; ~9 h for 114 species on hawkcpu).
#SBATCH --partition=CHANGE_ME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M
#SBATCH --time=12:00:00
#SBATCH --job-name=convgeno_throughput_calibration
#SBATCH --output=logs/throughput_calibration_%j.out
#SBATCH --error=logs/throughput_calibration_%j.err
#SBATCH --export=ALL

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "$PROJECT_ROOT"

# EDIT: activate the convgeno env so `diamond` and `python` are on PATH.
# The simplest, portable approach is to copy the conda bootstrap block from the
# top of your generated slurm_scripts/orthofinder_prepare.sh (between
# "set -euo pipefail" and the first echo) and paste it here -- it already has
# your cluster's exact module/conda paths. Placeholder:
#   module load <conda_module>            # if your cluster uses modules
#   source <conda_base>/etc/profile.d/conda.sh
#   conda activate <conda_env_prefix>

command -v diamond >/dev/null 2>&1 || {
    echo "ERROR: diamond not on PATH after env activation (edit the bootstrap)." >&2
    exit 1
}

CONFIG="$PROJECT_ROOT/pipeline_config.yaml"
CALIBRATOR="$PROJECT_ROOT/tools/throughput-constant/throughput_calibration.sh"

# Resolve the prepare outputs from the config.
OUT=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys
import yaml
print(yaml.safe_load(Path(sys.argv[1]).read_text())["orthofinder"]["output_dir"])
PY
)
PARENT=$(dirname "$OUT")
RUN=$(basename "$OUT")
SEARCH_COMMANDS="$PARENT/${RUN}_search_commands.txt"
WORK_DIR_FILE="$PARENT/${RUN}_working_dir_path.txt"
[ -f "$WORK_DIR_FILE" ] || { echo "ERROR: run prepare first ($WORK_DIR_FILE missing)." >&2; exit 1; }
WORK_DIR=$(<"$WORK_DIR_FILE")

echo "============================================================"
echo "DIAMOND throughput calibration"
echo "Job ID:           ${SLURM_JOB_ID:-?}"
echo "Node:             $HOSTNAME"
echo "Search commands:  $SEARCH_COMMANDS"
echo "WorkingDirectory: $WORK_DIR"
echo "Start:            $(date)"
echo "============================================================"

bash "$CALIBRATOR" "$SEARCH_COMMANDS" "$WORK_DIR" 100

echo "============================================================"
echo "Calibration finished. Copy the '>>> throughput_const = …' value into"
echo "pipeline_config.yaml under multinode.throughput_const."
echo "End: $(date)"
echo "============================================================"
