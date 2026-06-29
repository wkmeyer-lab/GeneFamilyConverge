"""
Tests for convgeno.slurm.runtime

Run with:  pytest tests/test_runtime.py -v
"""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from convgeno.slurm.runtime import (
    CondaRuntimeConfig,
    detect_conda_runtime,
    render_conda_bootstrap,
    runtime_config_from_dict,
    runtime_config_to_dict,
)


@pytest.fixture()
def sample_runtime() -> CondaRuntimeConfig:
    return CondaRuntimeConfig(
        conda_module="miniforge3/24.3.0-0",
        conda_base=Path("/share/apps/miniforge3/24.3.0-0"),
        conda_env_prefix=Path("/home/prm526/.conda/envs/convgeno"),
    )


@pytest.fixture()
def runtime_no_module() -> CondaRuntimeConfig:
    return CondaRuntimeConfig(
        conda_module=None,
        conda_base=Path("/opt/conda"),
        conda_env_prefix=Path("/home/user/envs/convgeno"),
    )


class TestRenderBootstrapContainsCondaSource:
    def test_uses_absolute_path_for_conda_sh(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert '[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]' in output
        assert 'source "$CONDA_BASE/etc/profile.d/conda.sh"' in output


class TestRenderBootstrapContainsActivatePrefix:
    def test_activates_by_absolute_prefix(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert 'conda activate "$CONDA_ENV"' in output


class TestRenderBootstrapWithModule:
    def test_includes_module_load(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert 'module load "$CONDA_MODULE"' in output

    def test_module_warns_on_failure(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert 'echo "WARNING: failed to load module: $CONDA_MODULE"' in output

    def test_module_no_silent_swallow(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "2>/dev/null || true" not in output


class TestRenderBootstrapWithoutModule:
    def test_no_module_load(self, runtime_no_module):
        output = render_conda_bootstrap(runtime_no_module)
        assert "module load" not in output

    def test_conda_module_empty_when_unset(self, runtime_no_module):
        output = render_conda_bootstrap(runtime_no_module)
        assert 'CONDA_MODULE=""' in output


class TestRenderBootstrapContainsDiagnostics:
    def test_which_python(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "which python" in output

    def test_which_orthofinder(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "which orthofinder" in output

    def test_which_diamond(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "which diamond" in output

    def test_which_convgeno(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "which convgeno" in output


class TestRenderBootstrapHasErrorChecks:
    def test_exit_127_for_conda_not_found(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "exit 127" in output

    def test_checks_conda_default_env(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "CONDA_DEFAULT_ENV" in output

    def test_exit_1_on_activation_failure(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "exit 1" in output


class TestRenderBootstrapNoCondaInfoBase:
    def test_conda_info_base_only_in_elif_fallback(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        primary_end = output.index("elif command -v conda")
        assert "conda info --base" not in output[:primary_end]

    def test_conda_info_base_in_elif_fallback(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "elif command -v conda >/dev/null 2>&1; then" in output
        assert 'source "$(conda info --base)/etc/profile.d/conda.sh"' in output

    def test_conda_info_base_only_in_elif_without_module(self, runtime_no_module):
        output = render_conda_bootstrap(runtime_no_module)
        primary_end = output.index("elif command -v conda")
        assert "conda info --base" not in output[:primary_end]


class TestRenderBootstrapGuards:
    def test_guards_conda_base(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert '[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]' in output

    def test_has_conda_base_fallback(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert "elif command -v conda >/dev/null 2>&1; then" in output

    def test_uses_bash_variables(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert 'CONDA_BASE="/share/apps/miniforge3/24.3.0-0"' in output
        assert 'CONDA_ENV="/home/prm526/.conda/envs/convgeno"' in output
        assert 'CONDA_MODULE="miniforge3/24.3.0-0"' in output
        assert "$CONDA_BASE" in output
        assert "$CONDA_ENV" in output


class TestRenderBootstrapMarkers:
    def test_starts_with_marker(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert output.startswith("# ---- convgeno runtime bootstrap ----")

    def test_ends_with_marker(self, sample_runtime):
        output = render_conda_bootstrap(sample_runtime)
        assert output.rstrip().endswith("# ---- end convgeno runtime bootstrap ----")


class TestRuntimeConfigRoundtrip:
    def test_serialize_deserialize_with_module(self, sample_runtime):
        d = runtime_config_to_dict(sample_runtime)
        restored = runtime_config_from_dict(d["runtime"])
        assert restored.conda_module == sample_runtime.conda_module
        assert restored.conda_base == sample_runtime.conda_base
        assert restored.conda_env_prefix == sample_runtime.conda_env_prefix

    def test_serialize_deserialize_without_module(self, runtime_no_module):
        d = runtime_config_to_dict(runtime_no_module)
        restored = runtime_config_from_dict(d["runtime"])
        assert restored.conda_module is None
        assert restored.conda_base == runtime_no_module.conda_base
        assert restored.conda_env_prefix == runtime_no_module.conda_env_prefix

    def test_dict_structure(self, sample_runtime):
        d = runtime_config_to_dict(sample_runtime)
        assert "runtime" in d
        assert d["runtime"]["conda_module"] == "miniforge3/24.3.0-0"
        assert d["runtime"]["conda_base"] == "/share/apps/miniforge3/24.3.0-0"
        assert d["runtime"]["conda_env_prefix"] == "/home/prm526/.conda/envs/convgeno"


class TestRuntimeConfigFromDictMissingKeys:
    def test_missing_conda_base_raises(self):
        with pytest.raises(ValueError, match="conda_base"):
            runtime_config_from_dict(
                {"conda_env_prefix": "/home/user/envs/convgeno"}
            )

    def test_missing_conda_env_prefix_raises(self):
        with pytest.raises(ValueError, match="conda_env_prefix"):
            runtime_config_from_dict({"conda_base": "/opt/conda"})

    def test_empty_dict_raises(self):
        with pytest.raises(ValueError):
            runtime_config_from_dict({})

    def test_conda_module_optional(self):
        config = runtime_config_from_dict({
            "conda_base": "/opt/conda",
            "conda_env_prefix": "/home/user/envs/convgeno",
        })
        assert config.conda_module is None


class TestDetectCondaRuntimeBasic:
    def test_detects_prefix_from_sys_prefix(self, tmp_path):
        fake_env = tmp_path / "envs" / "convgeno"
        fake_env.mkdir(parents=True)
        (fake_env / "conda-meta").mkdir()
        (fake_env / "bin").mkdir()
        (fake_env / "bin" / "python").touch()

        fake_base = tmp_path / "conda_base"
        fake_base.mkdir()
        (fake_base / "etc" / "profile.d").mkdir(parents=True)
        (fake_base / "etc" / "profile.d" / "conda.sh").touch()

        class FakeResult:
            returncode = 0
            stdout = str(fake_base) + "\n"
            stderr = ""

        with patch("sys.prefix", str(fake_env)), \
             patch("subprocess.run", return_value=FakeResult()), \
             patch.dict(os.environ, {"LOADEDMODULES": ""}, clear=False):
            config = detect_conda_runtime()

        assert config.conda_env_prefix == fake_env
        assert config.conda_base == fake_base
        assert config.conda_module is None

    def test_detects_conda_module_from_env(self, tmp_path):
        fake_env = tmp_path / "envs" / "convgeno"
        fake_env.mkdir(parents=True)
        (fake_env / "conda-meta").mkdir()
        (fake_env / "bin").mkdir()
        (fake_env / "bin" / "python").touch()

        fake_base = tmp_path / "conda_base"
        fake_base.mkdir()
        (fake_base / "etc" / "profile.d").mkdir(parents=True)
        (fake_base / "etc" / "profile.d" / "conda.sh").touch()

        class FakeResult:
            returncode = 0
            stdout = str(fake_base) + "\n"
            stderr = ""

        modules_str = "gcc/12.1.0:miniforge3/24.3.0-0:openmpi/4.1"

        with patch("sys.prefix", str(fake_env)), \
             patch("subprocess.run", return_value=FakeResult()), \
             patch.dict(os.environ, {"LOADEDMODULES": modules_str}, clear=False):
            config = detect_conda_runtime()

        assert config.conda_module == "miniforge3/24.3.0-0"

    def test_raises_when_conda_not_found(self, tmp_path):
        fake_env = tmp_path / "envs" / "convgeno"
        fake_env.mkdir(parents=True)
        (fake_env / "conda-meta").mkdir()

        with patch("sys.prefix", str(fake_env)), \
             patch("subprocess.run", side_effect=FileNotFoundError()):
            with pytest.raises(RuntimeError, match="conda command not found"):
                detect_conda_runtime()

    def test_raises_when_conda_info_fails(self, tmp_path):
        fake_env = tmp_path / "envs" / "convgeno"
        fake_env.mkdir(parents=True)
        (fake_env / "conda-meta").mkdir()

        class FakeResult:
            returncode = 1
            stdout = ""
            stderr = "some error"

        with patch("sys.prefix", str(fake_env)), \
             patch("subprocess.run", return_value=FakeResult()):
            with pytest.raises(RuntimeError, match="conda info --base failed"):
                detect_conda_runtime()
