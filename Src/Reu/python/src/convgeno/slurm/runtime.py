"""Centralized SLURM runtime environment bootstrap for compute nodes.

This module is the single source of truth for conda activation shell code.
No other module should contain inline conda bootstrap logic.

The core problem solved: Generated SLURM scripts must be fully self-contained
because compute nodes start with a bare shell — no conda on PATH, no modules
loaded, no user shell configuration sourced. Using ``$(conda info --base)``
in generated scripts is a bootstrap paradox (it requires conda to find conda).
Instead, we detect absolute paths at ``convgeno init`` time and embed them
directly into every generated script.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

__all__ = [
    "CondaRuntimeConfig",
    "detect_conda_runtime",
    "render_conda_bootstrap",
    "render_mafft_msa_shim",
    "render_ultrametric_autostep",
    "render_phenotype_tree_autostep",
    "runtime_config_to_dict",
    "runtime_config_from_dict",
]


def render_ultrametric_autostep(
    project_dir: str,
    results_expr: str,
    calibration_arg: str | None = None,
    *,
    user_tree: str | None = None,
    user_tree_is_ultrametric: bool = False,
    num_sites: int | None = None,
    species_dir: str | None = None,
) -> str:
    """Shell block that produces the ultrametric ``-t`` tree after a run.

    Emitted at the tail of a successful OrthoFinder job (single- or multi-node)
    so ``species_tree_ultrametric.nwk`` is produced automatically. It is
    **non-fatal** (OrthoFinder's own success is never undone by it) and, when it
    runs r8s, **skips cleanly** if r8s is not installed (``--skip-if-unavailable``).

    *results_expr* is the shell expression pointing at the OrthoFinder output or
    Results directory — ``$OUTPUT_DIR`` for single-node (results at
    ``$OUTPUT_DIR/Results_*``), ``$FINAL_RESULTS`` for multi-node (the deep
    ``$WORK_DIR/OrthoFinder/Results_*`` the resume script already resolves).
    ``make_tree_ultrametric.py`` handles either depth. It is only emitted when
    OrthoFinder output is actually needed (see below).

    *calibration_arg* is a ``make_tree_ultrametric.py --calibration`` value
    (``NAME:SP1,SP2:AGE``) from ``convgeno init``. When given, r8s scales the
    tree to real time; when ``None``, the root is anchored (relative time).

    User-supplied species tree (``convgeno init``'s ``species_tree`` block):

    - ``user_tree`` set + ``user_tree_is_ultrametric`` → the tree is taken as-is
      (``--assume-ultrametric``): validated and copied, **no r8s**, and no
      OrthoFinder output is read.
    - ``user_tree`` set, not ultrametric, ``num_sites`` given → r8s dates the
      user tree with that ``--nsites``; **no** OrthoFinder output is read.
    - ``user_tree`` set, not ultrametric, no ``num_sites`` → r8s dates the user
      tree but reads OrthoFinder's alignment (``results_expr``) for ``nsites``.
    - ``user_tree`` unset → r8s dates OrthoFinder's own ``SpeciesTree_rooted.txt``.

    *species_dir* (the proteome input dir) is passed as ``--species-dir`` for a
    user tree so its tips are re-checked against the species set at run time.
    """
    has_user_tree = bool(user_tree)
    # r8s runs unless a verified-ultrametric user tree was supplied.
    runs_r8s = not (has_user_tree and user_tree_is_ultrametric)
    # OrthoFinder output is needed to discover the tree (no user tree) or to
    # derive nsites (user tree without an explicit num_sites).
    needs_results = runs_r8s and (not has_user_tree or num_sites is None)

    args: list[str] = []
    if needs_results:
        args.append(f'"{results_expr}"')
    if has_user_tree:
        args.append(f'--input-tree "{user_tree}"')
        if user_tree_is_ultrametric:
            args.append("--assume-ultrametric")
        elif num_sites is not None:
            args.append(f"--nsites {num_sites}")
    args.append('-o "$OUTPUT_DIR/species_tree_ultrametric.nwk"')
    if runs_r8s:
        args.append(
            f'--calibration "{calibration_arg}"' if calibration_arg else "--root-age 1"
        )
    if has_user_tree and species_dir:
        args.append(f'--species-dir "{species_dir}"')
    if runs_r8s:
        args.append("--skip-if-unavailable")

    script = f"{project_dir}/Src/Loc/scripts/make_tree_ultrametric.py"
    command = " \\\n".join([f'python "{script}"'] + [f"    {a}" for a in args])

    if runs_r8s:
        setup = (
            'echo "Attempting automatic ultrametric conversion of the species '
            'tree (r8s)..."\n'
            'export PATH="$HOME/.local/bin:$PATH"   '
            "# r8s from tools/r8s/install_r8s.sh installs here"
        )
    else:
        setup = (
            'echo "Installing the user-supplied ultrametric species tree '
            '(validate + copy; no r8s)..."'
        )

    return f'''\
# ---- convgeno: auto ultrametric step (r8s; optional, non-fatal) ----
{setup}
set +e
{command}
CONVGENO_R8S_STATUS=$?
set -e
if [ "$CONVGENO_R8S_STATUS" -ne 0 ]; then
    echo "WARNING: ultrametric step returned $CONVGENO_R8S_STATUS (non-fatal)."
fi'''


def render_phenotype_tree_autostep(
    project_dir: str,
    *,
    phenotype_table: str,
    id_col: str = "species",
    pheno_col: str = "phenotype",
    model: str = "ER",
) -> str:
    """Shell block that builds the categorical phenotype tree (CAFE ``-y``).

    Emitted at the tail of a successful OrthoFinder job, immediately AFTER the
    ultrametric autostep, so the CAFE ``-y`` multi-lambda tree is produced
    automatically as soon as the ultrametric ``-t`` tree and the phenotype tip
    data exist. It reconstructs ancestral phenotype states on
    ``$OUTPUT_DIR/species_tree_ultrametric.nwk`` (the r8s output) and writes
    ``$OUTPUT_DIR/lambda_tree.nwk`` (+ ``lambda_legend.tsv``) beside it.

    Like the r8s step it is **non-fatal** (OrthoFinder's success is never undone)
    and **skips cleanly** when it cannot run: no ``Rscript`` on PATH (R + the
    ``tools/r-deps`` packages absent), no ultrametric tree yet, or the phenotype
    table missing. It is only emitted when a phenotype table was configured at
    ``convgeno init`` (``PhenotypeTreeConfig``). The script is invoked with
    explicit flags — no config file is read at run time — so the generated SLURM
    script stays fully self-contained.
    """
    script = f"{project_dir}/Src/Loc/scripts/make_categorical_phenotype_tree.R"
    args = [
        '--tree "$OUTPUT_DIR/species_tree_ultrametric.nwk"',
        f'--phenotypes "{phenotype_table}"',
        f'--id-col "{id_col}"',
        f'--pheno-col "{pheno_col}"',
        f'--model "{model}"',
        '--out "$OUTPUT_DIR/lambda_tree.nwk"',
        '--legend "$OUTPUT_DIR/lambda_legend.tsv"',
        f'--reu-dir "{project_dir}/Src/Reu/r"',
    ]
    command = " \\\n".join([f'Rscript "{script}"'] + [f"    {a}" for a in args])

    return f'''\
# ---- convgeno: auto categorical phenotype tree (CAFE -y; optional, non-fatal) ----
echo "Attempting automatic categorical phenotype tree (CAFE -y)..."
if ! command -v Rscript >/dev/null 2>&1; then
    echo "WARNING: Rscript not found; skipping phenotype tree (need R + tools/r-deps)."
elif [ ! -s "$OUTPUT_DIR/species_tree_ultrametric.nwk" ]; then
    echo "WARNING: no ultrametric tree; skipping phenotype tree."
elif [ ! -f "{phenotype_table}" ]; then
    echo "WARNING: phenotype table not found; skipping phenotype tree."
else
    set +e
    {command}
    CONVGENO_PHENO_STATUS=$?
    set -e
    if [ "$CONVGENO_PHENO_STATUS" -ne 0 ]; then
        echo "WARNING: phenotype tree step exit $CONVGENO_PHENO_STATUS (non-fatal)."
    fi
fi'''


def render_mafft_msa_shim() -> str:
    """Shell block that shadows ``mafft`` on PATH with a fast, retrying shim.

    Root cause it addresses: OrthoFinder's default MSA command for orthogroups
    with < 500 sequences is MAFFT **L-INS-i** (``mafft --localpair --maxiterate
    1000 --anysymbol``), which is pathologically slow on moderately large gene
    families and, under the multi-node resume, sporadically produces an **empty**
    alignment. An empty alignment for a single-copy orthogroup is fatal
    ("Species tree inference failed"). The fast command (``mafft --anysymbol``,
    which OrthoFinder itself already uses for >= 500-sequence OGs) aligns the same
    sequences instantly and reliably.

    This shim forces the fast command for **every** orthogroup, with up to 3
    retries and a non-empty-output check, so a stray failure can never leave an
    empty alignment. Non-alignment invocations (the startup dependency test,
    ``--version``) pass straight through to the real mafft.

    It is installed into a fresh ``mktemp -d`` directory prepended to ``PATH`` --
    job-local, no ``$HOME``/config global state, and no edits to OrthoFinder's
    ``config.json``. Emit it AFTER the conda bootstrap (so the real ``mafft`` is
    resolvable) and BEFORE invoking ``orthofinder``.
    """
    return '''\
# ---- convgeno MAFFT MSA shim: fast --anysymbol + retry (no L-INS-i) ----
CONVGENO_REAL_MAFFT="$(command -v mafft)"
export CONVGENO_REAL_MAFFT
CONVGENO_SHIM_DIR="$(mktemp -d)"
cat > "$CONVGENO_SHIM_DIR/mafft" <<'CONVGENO_MAFFT_SHIM'
#!/bin/bash
# convgeno MAFFT MSA shim (auto-generated; see runtime.render_mafft_msa_shim).
set -u
REAL="${CONVGENO_REAL_MAFFT:-mafft}"
INPUT="${@: -1}"
# Non-alignment call (dependency test / --version / no input file): pass through.
if [ "$#" -eq 0 ] || [ ! -f "$INPUT" ]; then
    exec "$REAL" "$@"
fi
# Alignment call: force the fast command, retry, require non-empty output.
tmp="$(mktemp)"
for _attempt in 1 2 3; do
    if "$REAL" --anysymbol "$INPUT" > "$tmp" 2>/dev/null && [ -s "$tmp" ]; then
        cat "$tmp"
        rm -f "$tmp"
        exit 0
    fi
done
rm -f "$tmp"
echo "convgeno mafft shim: alignment failed for $INPUT after 3 attempts" >&2
exit 1
CONVGENO_MAFFT_SHIM
chmod +x "$CONVGENO_SHIM_DIR/mafft"
export PATH="$CONVGENO_SHIM_DIR:$PATH"
echo "MAFFT shim active (fast --anysymbol + retry); real: $CONVGENO_REAL_MAFFT"'''


@dataclass(frozen=True)
class CondaRuntimeConfig:
    """Absolute paths needed to bootstrap conda on a bare compute node.

    Attributes
    ----------
    conda_module
        HPC module to load before conda is available (e.g.
        ``"miniforge3/24.3.0-0"``). ``None`` if the cluster does not
        use environment modules or conda is available without one.
    conda_base
        Absolute path to the conda installation directory. Used to
        source ``etc/profile.d/conda.sh`` without requiring conda
        on PATH first.
    conda_env_prefix
        Absolute path to the target conda environment. More robust
        than activating by name because it does not depend on
        ``~/.condarc`` search paths.
    """

    conda_module: Optional[str]
    conda_base: Path
    conda_env_prefix: Path


def detect_conda_runtime() -> CondaRuntimeConfig:
    """Detect conda runtime configuration from the current environment.

    Must be called while the ``convgeno`` conda environment is active
    (i.e. during ``convgeno init`` on a login node where conda works).

    Returns
    -------
    CondaRuntimeConfig
        Populated configuration ready for serialization.

    Raises
    ------
    RuntimeError
        If ``conda info --base`` fails (conda not available in the
        current shell).
    """
    conda_env_prefix = Path(sys.prefix)

    conda_meta = conda_env_prefix / "conda-meta"
    if not conda_meta.is_dir():
        logger.warning(
            "WARNING: convgeno init does not appear to be running inside "
            "a conda environment."
        )

    try:
        result = subprocess.run(
            ["conda", "info", "--base"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode != 0:
            raise RuntimeError(
                "conda info --base failed. Is conda available in your shell?\n"
                f"stderr: {result.stderr.strip()}"
            )
        conda_base = Path(result.stdout.strip())
    except FileNotFoundError:
        raise RuntimeError(
            "conda command not found. Activate your conda environment "
            "before running 'convgeno init'."
        )

    conda_sh = conda_base / "etc" / "profile.d" / "conda.sh"
    if not conda_sh.exists():
        logger.warning(
            "WARNING: conda.sh not found at expected path: %s", conda_sh
        )

    python_bin = conda_env_prefix / "bin" / "python"
    if not python_bin.exists():
        python_bin_alt = conda_env_prefix / "python.exe"
        if not python_bin_alt.exists():
            logger.warning(
                "WARNING: python not found in environment: %s",
                conda_env_prefix,
            )

    orthofinder_bin = conda_env_prefix / "bin" / "orthofinder"
    if not orthofinder_bin.exists():
        logger.warning(
            "WARNING: orthofinder not found in environment %s. "
            "Install it before running the pipeline.",
            conda_env_prefix,
        )

    conda_module = _detect_conda_module()

    return CondaRuntimeConfig(
        conda_module=conda_module,
        conda_base=conda_base,
        conda_env_prefix=conda_env_prefix,
    )


def _detect_conda_module() -> Optional[str]:
    """Check LOADEDMODULES for a conda-related module name.

    Looks for module names containing 'conda', 'miniforge',
    'miniconda', or 'anaconda'.
    """
    loaded = os.environ.get("LOADEDMODULES", "")
    if not loaded:
        return None

    keywords = ("conda", "miniforge", "miniconda", "anaconda")
    for module_name in loaded.split(":"):
        lower = module_name.lower()
        if any(kw in lower for kw in keywords):
            return module_name

    return None


def render_conda_bootstrap(runtime: CondaRuntimeConfig) -> str:
    """Render the bash bootstrap block for a SLURM script.

    The returned string should be inserted after ``set -euo pipefail``
    and before any pipeline commands. It does NOT include
    ``set -euo pipefail`` or ``#SBATCH`` directives — those are the
    caller's responsibility.

    Parameters
    ----------
    runtime
        The conda runtime configuration (absolute paths, optional
        module name).

    Returns
    -------
    str
        Multi-line bash code that bootstraps the environment on a
        bare compute node.
    """
    conda_base_str = runtime.conda_base.as_posix()
    env_prefix_str = runtime.conda_env_prefix.as_posix()
    conda_module_str = runtime.conda_module or ""
    env_name = runtime.conda_env_prefix.name

    lines: list[str] = []
    lines.append("# ---- convgeno runtime bootstrap ----")
    lines.append("")
    lines.append('echo "convgeno runtime bootstrap"')
    lines.append('echo "  Date:     $(date)"')
    lines.append('echo "  Hostname: $(hostname)"')
    lines.append('echo "  Job ID:   ${SLURM_JOB_ID:-none}"')
    lines.append("")
    lines.append(f'CONDA_MODULE="{conda_module_str}"')
    lines.append(f'CONDA_BASE="{conda_base_str}"')
    lines.append(f'CONDA_ENV="{env_prefix_str}"')
    lines.append("")

    if runtime.conda_module is not None:
        lines.append("# Load conda module if configured")
        lines.append('if [ -n "$CONDA_MODULE" ]; then')
        lines.append("    if command -v module >/dev/null 2>&1; then")
        lines.append('        module load "$CONDA_MODULE" || {')
        lines.append(
            '            echo "WARNING: failed to load module: $CONDA_MODULE" >&2'
        )
        lines.append("        }")
        lines.append("    fi")
        lines.append("fi")
        lines.append("")

    lines.append("# Source conda shell hook using absolute path")
    lines.append('if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then')
    lines.append('    source "$CONDA_BASE/etc/profile.d/conda.sh"')
    lines.append("elif command -v conda >/dev/null 2>&1; then")
    lines.append('    source "$(conda info --base)/etc/profile.d/conda.sh"')
    lines.append("else")
    lines.append(
        '    echo "ERROR: could not initialize conda. '
        'Check runtime.conda_base in pipeline_config.yaml." >&2'
    )
    lines.append("    exit 127")
    lines.append("fi")
    lines.append("")

    lines.append("# Verify conda is now available")
    lines.append("if ! command -v conda &> /dev/null; then")
    lines.append('    echo "ERROR: conda command not available after sourcing conda.sh"')
    lines.append('    echo "Check that the conda installation at $CONDA_BASE is intact."')
    lines.append("    exit 127")
    lines.append("fi")
    lines.append("")

    lines.append("# Activate environment by absolute prefix")
    lines.append('conda activate "$CONDA_ENV"')
    lines.append("")

    lines.append("# Verify activation succeeded")
    lines.append(f'if [ "$CONDA_DEFAULT_ENV" != "{env_name}" ] && '
                 f'[ "$CONDA_PREFIX" != "$CONDA_ENV" ]; then')
    lines.append('    echo "ERROR: Failed to activate conda environment at $CONDA_ENV"')
    lines.append('    echo "  CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"')
    lines.append('    echo "  CONDA_PREFIX=$CONDA_PREFIX"')
    lines.append("    exit 1")
    lines.append("fi")
    lines.append("")

    lines.append("# Diagnostic info")
    lines.append('echo "  Python:      $(which python || true)"')
    lines.append('echo "  OrthoFinder: $(which orthofinder || true)"')
    lines.append('echo "  DIAMOND:     $(which diamond || true)"')
    lines.append('echo "  convgeno:    $(which convgeno || true)"')
    lines.append('echo ""')
    lines.append("")
    lines.append("# ---- end convgeno runtime bootstrap ----")

    return "\n".join(lines)


def runtime_config_to_dict(config: CondaRuntimeConfig) -> dict:
    """Serialize a CondaRuntimeConfig to a dict suitable for YAML output.

    Returns
    -------
    dict
        Structure: ``{"runtime": {"conda_module": ..., "conda_base": ...,
        "conda_env_prefix": ...}}``.
    """
    return {
        "runtime": {
            "conda_module": config.conda_module,
            "conda_base": config.conda_base.as_posix(),
            "conda_env_prefix": config.conda_env_prefix.as_posix(),
        }
    }


def runtime_config_from_dict(data: dict) -> CondaRuntimeConfig:
    """Deserialize a CondaRuntimeConfig from a YAML-loaded dict.

    Parameters
    ----------
    data
        The ``runtime`` section from the config file (the value under
        the ``"runtime"`` key, not the top-level dict).

    Raises
    ------
    ValueError
        If required keys are missing.
    """
    if "conda_base" not in data:
        raise ValueError(
            "Runtime config is missing required key 'conda_base'. "
            "Re-run 'convgeno init' to regenerate the config."
        )
    if "conda_env_prefix" not in data:
        raise ValueError(
            "Runtime config is missing required key 'conda_env_prefix'. "
            "Re-run 'convgeno init' to regenerate the config."
        )

    return CondaRuntimeConfig(
        conda_module=data.get("conda_module"),
        conda_base=Path(data["conda_base"]),
        conda_env_prefix=Path(data["conda_env_prefix"]),
    )
