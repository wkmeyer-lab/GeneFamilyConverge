"""
Tests for convgeno.slurm.script_generator

Run with:  pytest tests/test_script_generator.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig
from convgeno.slurm.script_generator import generate_orthofinder_script, write_script


@pytest.fixture()
def sample_config() -> PipelineConfig:
    return PipelineConfig(
        project_dir="/share/ceph/project",
        conda_env="convgeno",
        slurm=SlurmConfig(
            partition="hawkcpu", cpus_per_task=16, time_limit="48:00:00"
        ),
        orthofinder=OrthoFinderConfig(
            input_dir="/share/ceph/project/Data/interim/cleaned_proteomes",
            output_dir="/share/ceph/project/Data/processed/orthofinder",
        ),
    )


class TestGenerateOrthoFinderScript:
    def test_starts_with_shebang(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert script.startswith("#!/bin/bash\n")

    def test_contains_sbatch_directives(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "#SBATCH --partition=hawkcpu" in script
        assert "#SBATCH --time=48:00:00" in script
        assert "#SBATCH --cpus-per-task=16" in script
        assert "#SBATCH --job-name=convgeno_orthofinder" in script

    def test_contains_conda_activation(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "conda activate convgeno" in script

    def test_contains_orthofinder_command(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "orthofinder" in script
        assert "-f /share/ceph/project/Data/interim/cleaned_proteomes" in script
        assert "-o /share/ceph/project/Data/processed/orthofinder" in script
        assert "-t 16" in script
        assert "-a 8" in script
        assert "-S diamond" in script

    def test_contains_input_validation(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert 'if [ ! -d "$INPUT_DIR" ]' in script
        assert "FASTA_COUNT" in script

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

    def test_raises_without_orthofinder_config(self):
        config = PipelineConfig(
            project_dir="/project", slurm=SlurmConfig(partition="hawkcpu")
        )
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_orthofinder_script(config)

    def test_with_extra_orthofinder_args(self):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(
                input_dir="/in",
                output_dir="/out",
                extra_args=["--fewer-files"],
            ),
        )
        script = generate_orthofinder_script(config)
        assert "--fewer-files" in script


class TestWriteScript:
    def test_creates_file(self, tmp_path: Path):
        script_path = tmp_path / "test_job.sh"
        write_script("#!/bin/bash\necho hello", script_path)
        assert script_path.exists()
        assert script_path.read_text(encoding="utf-8") == "#!/bin/bash\necho hello"

    def test_is_executable(self, tmp_path: Path):
        script_path = tmp_path / "job.sh"
        write_script("#!/bin/bash\necho test", script_path)
        assert script_path.stat().st_mode & 0o111

    def test_creates_parent_dirs(self, tmp_path: Path):
        script_path = tmp_path / "subdir" / "nested" / "job.sh"
        write_script("#!/bin/bash", script_path)
        assert script_path.exists()
