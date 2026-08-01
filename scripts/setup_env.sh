#!/usr/bin/env bash
# Create or update the `convgeno` conda environment safely on an HPC cluster.
#
# WHY THIS SCRIPT EXISTS
#   On a login node, conda's dependency solve for this environment is often
#   OOM-killed by the node's memory cgroup. The symptom is:
#
#       Collecting package metadata (repodata.json): - Killed
#
#   That is NOT a dependency conflict (a real conflict prints a solver error, not
#   "Killed"); it is the login node terminating a memory-heavy solve. The fix is
#   to run the solve on a compute node (via `srun`) and use the lower-memory
#   libmamba solver. This script does both automatically: run it from a login
#   node and it grabs a compute node for you, then creates (or updates) the env.
#
# USAGE (from the repository root, or anywhere -- paths are resolved):
#   ./scripts/setup_env.sh                       # create, or update if it exists
#   ./scripts/setup_env.sh --here                # do NOT srun; solve on this node
#                                                #   (use inside an interactive
#                                                #    allocation, or on a laptop)
#   ./scripts/setup_env.sh --mem 32G --time 03:00:00 --cpus 4 --partition <p>
#   ./scripts/setup_env.sh --name myenv --conda-module miniforge3/24.3.0-0
#
# ENV VAR OVERRIDES: ENV_NAME, SETUP_MEM, SETUP_TIME, SETUP_CPUS,
#                    SETUP_PARTITION, CONDA_MODULE.
#
# After it finishes, activate and run the pipeline:
#   conda activate convgeno
#   convgeno init && convgeno run

set -Eeuo pipefail

# --- resolve absolute paths (survive the srun re-exec, any cwd) --------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_FILE="${REPO_ROOT}/environment.yml"

# --- defaults (overridable by env var, then by flag) -------------------------
ENV_NAME="${ENV_NAME:-convgeno}"
MEM="${SETUP_MEM:-32G}"
TIME="${SETUP_TIME:-02:00:00}"
CPUS="${SETUP_CPUS:-4}"
PARTITION="${SETUP_PARTITION:-}"
CONDA_MODULE="${CONDA_MODULE:-miniconda3}"
HERE=false

log()  { printf '[setup_env] %s\n' "$*"; }
warn() { printf '[setup_env] WARNING: %s\n' "$*" >&2; }
fail() { printf '[setup_env] ERROR: %s\n' "$*" >&2; exit 1; }

# --- parse flags -------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --here)          HERE=true; shift ;;
        --name)          ENV_NAME="$2"; shift 2 ;;
        --mem)           MEM="$2"; shift 2 ;;
        --time)          TIME="$2"; shift 2 ;;
        --cpus)          CPUS="$2"; shift 2 ;;
        --partition)     PARTITION="$2"; shift 2 ;;
        --conda-module)  CONDA_MODULE="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,30p' "$SCRIPT_PATH"; exit 0 ;;
        *) fail "Unknown argument: $1 (try --help)" ;;
    esac
done

[[ -f "$ENV_FILE" ]] || fail "environment.yml not found at: $ENV_FILE"

# --- decide whether to hop to a compute node ---------------------------------
# Skip the hop if: the user forced --here, we're already inside a SLURM
# allocation, or there's no srun at all (local/dev machine).
need_srun() {
    $HERE && return 1
    [[ -n "${SLURM_JOB_ID:-}" ]] && return 1
    command -v srun >/dev/null 2>&1 || return 1
    return 0
}

if need_srun; then
    SRUN_ARGS=(
        --job-name=convgeno_env_setup
        --cpus-per-task="$CPUS"
        --mem="$MEM"
        --time="$TIME"
    )
    [[ -n "$PARTITION" ]] && SRUN_ARGS+=(--partition="$PARTITION")
    # Only request a pseudo-terminal when we actually have one, so this still
    # works when stdout is redirected or run non-interactively.
    [[ -t 1 ]] && SRUN_ARGS+=(--pty)
    log "Login/submit node detected. Moving the conda solve onto a compute node"
    log "  (avoids the login-node OOM kill). Requesting: cpus=$CPUS mem=$MEM time=$TIME${PARTITION:+ partition=$PARTITION}"
    log "  srun ${SRUN_ARGS[*]} bash -l $SCRIPT_PATH --here"
    # bash -l so the compute-node login profile is sourced and `module` works.
    # --here prevents a second hop inside the allocation.
    exec srun "${SRUN_ARGS[@]}" bash -l "$SCRIPT_PATH" --here \
        --name "$ENV_NAME" --conda-module "$CONDA_MODULE"
fi

# ----------------------------------------------------------------------------
# From here on we are on a node where it's safe to solve.
# ----------------------------------------------------------------------------
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    log "Running inside SLURM allocation ${SLURM_JOB_ID} on $(hostname)."
else
    log "Running directly on $(hostname) (no SLURM hop)."
fi

# --- make conda available ----------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
    log "conda not on PATH; attempting: module load $CONDA_MODULE"
    if command -v module >/dev/null 2>&1; then
        module load "$CONDA_MODULE" 2>/dev/null \
            || warn "Could not 'module load $CONDA_MODULE'. Override with --conda-module."
    else
        warn "No 'module' command available to load conda."
    fi
fi
command -v conda >/dev/null 2>&1 \
    || fail "conda is not available. Load your conda/miniconda module first, or pass --conda-module."

log "Using conda: $(command -v conda)  ($(conda --version 2>/dev/null || echo 'unknown version'))"

# --- prefer the low-memory libmamba solver if this conda supports it ---------
SOLVER_ARGS=()
if conda env create --help 2>/dev/null | grep -q -- '--solver'; then
    SOLVER_ARGS=(--solver libmamba)
    log "Using the libmamba solver (lower memory than the classic solver)."
else
    warn "This conda has no --solver flag; using the default solver. If the solve"
    warn "is killed, upgrade conda or 'conda install -n base conda-libmamba-solver'."
fi

# --- create (fresh) or update (existing), from the repo root -----------------
# CWD must be the repo root: environment.yml pins the package with an editable
# pip path ('-e ./Src/Reu/python/') that conda resolves relative to CWD.
cd "$REPO_ROOT"

env_exists() { conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "$ENV_NAME"; }

if env_exists; then
    log "Environment '$ENV_NAME' already exists -> updating (conda env update --prune)."
    conda env update -n "$ENV_NAME" -f "$ENV_FILE" --prune "${SOLVER_ARGS[@]}"
else
    log "Creating environment '$ENV_NAME' from $ENV_FILE."
    conda env create -n "$ENV_NAME" -f "$ENV_FILE" "${SOLVER_ARGS[@]}"
fi

log "conda solve finished successfully."

# --- verify (non-fatal) ------------------------------------------------------
log "Verifying the environment ..."
CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    if conda activate "$ENV_NAME" 2>/dev/null; then
        printf '  python:     %s\n' "$(python --version 2>&1 || echo 'MISSING')"
        if command -v orthofinder >/dev/null 2>&1; then printf '  orthofinder: OK\n'; else printf '  orthofinder: MISSING\n'; fi
        printf '  snakemake:  %s\n' "$(snakemake --version 2>/dev/null || echo 'MISSING')"
        Rscript -e 'ok<-TRUE; for (p in c("ape","castor","expm","yaml")) if(!requireNamespace(p,quietly=TRUE)){cat("  R pkg MISSING:",p,"\n"); ok<-FALSE}; if(ok) cat("  R deps:     ape/castor/expm/yaml OK\n")' 2>/dev/null \
            || printf '  R deps:     could not verify (Rscript failed)\n'
        # r8s and CAFE-5 are external tools (not on conda); just report presence.
        if command -v r8s   >/dev/null 2>&1; then printf '  r8s:        on PATH\n';   else printf '  r8s:        not on PATH (build with tools/r8s/install_r8s.sh; the dating step will skip)\n'; fi
        if command -v cafe5 >/dev/null 2>&1; then printf '  cafe5:      on PATH\n';    else printf '  cafe5:      not on PATH (install CAFE-5 separately; required for the turnover step)\n'; fi
    else
        warn "Could not activate '$ENV_NAME' for verification (the env was still built)."
    fi
else
    warn "conda profile script not found; skipping verification."
fi

cat <<EOF

[setup_env] Done. Next:

    conda activate ${ENV_NAME}
    convgeno init      # one-time interactive setup
    convgeno run       # run the full pipeline as a DAG

EOF
