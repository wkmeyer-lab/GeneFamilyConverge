"""
Tests for convgeno.slurm.script_generator

Run with:  pytest tests/test_script_generator.py -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig
from convgeno.slurm.runtime import CondaRuntimeConfig
from convgeno.slurm.script_generator import generate_orthofinder_script, write_script


@pytest.fixture()
def sample_runtime() -> CondaRuntimeConfig:
    return CondaRuntimeConfig(
        conda_module="miniforge3/24.3.0-0",
        conda_base=Path("/share/apps/miniforge3/24.3.0-0"),
        conda_env_prefix=Path("/home/prm526/.conda/envs/convgeno"),
    )


@pytest.fixture()
def sample_config(sample_runtime) -> PipelineConfig:
    return PipelineConfig(
        project_dir="/share/ceph/project",
        conda_env="convgeno",
        slurm=SlurmConfig(
            partition="hawkcpu", cpus_per_task=16, time_limit="72:00:00"
        ),
        orthofinder=OrthoFinderConfig(
            input_dir="/share/ceph/project/Data/interim/cleaned_proteomes",
            output_dir="/share/ceph/project/Data/processed/orthofinder",
        ),
        runtime=sample_runtime,
    )


class TestGenerateOrthoFinderScript:
    def test_starts_with_shebang(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert script.startswith("#!/bin/bash\n")

    def test_contains_sbatch_directives(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "#SBATCH --partition=hawkcpu" in script
        assert "#SBATCH --time=72:00:00" in script
        assert "#SBATCH --cpus-per-task=16" in script
        assert "#SBATCH --job-name=convgeno_orthofinder" in script

    def test_contains_conda_activation(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert 'conda activate "$CONDA_ENV"' in script

    def test_contains_orthofinder_command(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "orthofinder" in script
        assert "-f /share/ceph/project/Data/interim/cleaned_proteomes" in script
        assert "-o /share/ceph/project/Data/processed/orthofinder" in script
        assert "-t 16" in script
        assert "-a 8" in script
        assert "-S diamond" in script
        # OrthoFinder's -M takes a method, not a program name.
        assert "-M msa" in script
        assert "-A mafft" in script
        assert "-T fasttree" in script
        assert "-M mafft" not in script

    def test_contains_input_validation(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert 'if [ ! -d "$INPUT_DIR" ]' in script
        assert "FASTA_COUNT" in script

    def test_does_not_create_output_dir_directly(self, sample_config: PipelineConfig):
        # OrthoFinder refuses to run with an existing -o output directory.
        # The script must create only the PARENT directory.
        script = generate_orthofinder_script(sample_config)
        assert 'mkdir -p "$(dirname "$OUTPUT_DIR")"' in script
        assert 'mkdir -p "$OUTPUT_DIR"' not in script

    def test_fails_loudly_if_output_dir_exists(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert 'if [ -e "$OUTPUT_DIR" ]' in script
        assert "already exists" in script

    def test_uses_help_for_version_probe_not_unsupported_flag(
        self, sample_config: PipelineConfig
    ):
        # OrthoFinder 2.5.5 does not support --version. Use -h instead.
        script = generate_orthofinder_script(sample_config)
        assert "orthofinder --version" not in script
        assert "orthofinder -h" in script

    def test_contains_exit_code_check(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "EXIT_CODE=$?" in script
        assert "exit $EXIT_CODE" in script

    def test_contains_timing(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "START_SECONDS" in script
        assert "Duration:" in script

    def test_set_euo_pipefail(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "set -euo pipefail" in script

    def test_raises_without_orthofinder_config(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            runtime=sample_runtime,
        )
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_orthofinder_script(config)

    def test_with_extra_orthofinder_args(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(
                input_dir="/in",
                output_dir="/out",
                extra_args=["--fewer-files"],
            ),
            runtime=sample_runtime,
        )
        script = generate_orthofinder_script(config)
        assert "--fewer-files" in script

    def test_generated_script_contains_export_all(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        assert "#SBATCH --export=ALL" in script

    def test_generated_script_no_conda_info_base_as_primary(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        primary_end = script.index("elif command -v conda")
        assert "conda info --base" not in script[:primary_end]

    def test_contains_bootstrap_markers(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        assert "# ---- convgeno runtime bootstrap ----" in script
        assert "# ---- end convgeno runtime bootstrap ----" in script

    def test_sources_conda_sh_with_absolute_path(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        assert '[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]' in script
        assert 'source "$CONDA_BASE/etc/profile.d/conda.sh"' in script

    def test_raises_without_runtime(self):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=None,
        )
        with pytest.raises(ValueError, match="runtime"):
            generate_orthofinder_script(config)

    def test_explicit_runtime_overrides_config(self):
        alt_runtime = CondaRuntimeConfig(
            conda_module=None,
            conda_base=Path("/alt/conda"),
            conda_env_prefix=Path("/alt/envs/myenv"),
        )
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=CondaRuntimeConfig(
                conda_module="old/1.0",
                conda_base=Path("/old/conda"),
                conda_env_prefix=Path("/old/envs/convgeno"),
            ),
        )
        script = generate_orthofinder_script(config, alt_runtime)
        assert 'CONDA_BASE="/alt/conda"' in script
        assert 'conda activate "$CONDA_ENV"' in script


class TestWriteScript:
    def test_creates_file(self, tmp_path: Path):
        script_path = tmp_path / "test_job.sh"
        write_script("#!/bin/bash\necho hello", script_path)
        assert script_path.exists()
        assert script_path.read_text(encoding="utf-8") == "#!/bin/bash\necho hello"

    @pytest.mark.skipif(
        sys.platform == "win32",
        reason="Windows does not support Unix file permission bits",
    )
    def test_is_executable(self, tmp_path: Path):
        script_path = tmp_path / "job.sh"
        write_script("#!/bin/bash\necho test", script_path)
        assert script_path.stat().st_mode & 0o111

    def test_creates_parent_dirs(self, tmp_path: Path):
        script_path = tmp_path / "subdir" / "nested" / "job.sh"
        write_script("#!/bin/bash", script_path)
        assert script_path.exists()
