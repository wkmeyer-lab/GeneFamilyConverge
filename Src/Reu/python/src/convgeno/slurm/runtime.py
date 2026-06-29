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
    "runtime_config_to_dict",
    "runtime_config_from_dict",
]


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
