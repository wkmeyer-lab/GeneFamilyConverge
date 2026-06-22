"""
Tests for convgeno.slurm.multinode_generator

Run with:  pytest tests/test_multinode_generator.py -v
"""

from __future__ import annotations

import pytest

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig
from convgeno.slurm.multinode_generator import (
    generate_prepare_script,
    generate_resume_script,
    generate_search_array_script,
)


@pytest.fixture()
def sample_config() -> PipelineConfig:
    return PipelineConfig(
        project_dir="/share/ceph/project",
        conda_env="convgeno",
        slurm=SlurmConfig(
            partition="hawkcpu",
            cpus_per_task=16,
            time_limit="48:00:00",
            account="wym219",
        ),
        orthofinder=OrthoFinderConfig(
            input_dir="/share/ceph/project/proteomes",
            output_dir="/share/ceph/project/results",
            search_threads=16,
            analysis_threads=8,
        ),
    )


class TestPrepareScript:
    def test_shebang(self, sample_config):
        assert generate_prepare_script(sample_config).startswith("#!/bin/bash\n")

    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_prepare" in generate_prepare_script(sample_config)

    def test_single_node(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert "--nodes=1" in script
        assert "--cpus-per-task=4" in script

    def test_short_time_limit(self, sample_config):
        assert "--time=02:00:00" in generate_prepare_script(sample_config)

    def test_runs_op_flag(self, sample_config):
        assert "-op" in generate_prepare_script(sample_config)

    def test_captures_commands(self, sample_config):
        assert "diamond_commands.txt" in generate_prepare_script(sample_config)

    def test_uses_partition_from_config(self, sample_config):
        assert "--partition=hawkcpu" in generate_prepare_script(sample_config)

    def test_includes_account(self, sample_config):
        assert "--account=wym219" in generate_prepare_script(sample_config)

    def test_does_not_create_output_dir_directly(self, sample_config):
        # OrthoFinder refuses to run with an existing -o directory.
        # Only the parent should be created.
        script = generate_prepare_script(sample_config)
        assert 'mkdir -p "$(dirname "$OUTPUT_DIR")"' in script
        assert 'mkdir -p "$OUTPUT_DIR"' not in script

    def test_fails_loudly_if_output_dir_exists(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert 'if [ -e "$OUTPUT_DIR" ]' in script
        assert "already exists" in script


def _search_script(config, **kwargs):
    """Helper: call generate_search_array_script with sensible defaults."""
    kwargs.setdefault("commands_per_task", 50)
    kwargs.setdefault("array_max", 99)
    return generate_search_array_script(config, **kwargs)


class TestSearchArrayScript:
    def test_is_array_job(self, sample_config):
        assert "--array=" in _search_script(sample_config)

    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_search" in _search_script(sample_config)

    def test_reads_commands_file(self, sample_config):
        assert "diamond_commands.txt" in _search_script(sample_config)

    def test_uses_slurm_array_task_id(self, sample_config):
        assert "SLURM_ARRAY_TASK_ID" in _search_script(sample_config)

    def test_exits_if_no_work(self, sample_config):
        assert "exit 0" in _search_script(sample_config)

    def test_counts_failures(self, sample_config):
        assert "FAILED" in _search_script(sample_config)

    def test_uses_search_threads(self, sample_config):
        assert "--cpus-per-task=16" in _search_script(sample_config)

    def test_custom_commands_per_task(self, sample_config):
        script = _search_script(sample_config, commands_per_task=100)
        assert "COMMANDS_PER_TASK=100" in script

    def test_uses_computed_array_max(self, sample_config):
        # The hardcoded 0-9999 was rejected by SLURM clusters with
        # "Invalid job array specification". The array range must come
        # from caller-supplied array_max.
        script = generate_search_array_script(
            sample_config, commands_per_task=50, array_max=259
        )
        assert "#SBATCH --array=0-259" in script
        assert "#SBATCH --array=0-9999" not in script

    def test_array_max_zero_is_valid_single_task(self, sample_config):
        script = generate_search_array_script(
            sample_config, commands_per_task=50, array_max=0
        )
        assert "#SBATCH --array=0-0" in script

    def test_invalid_commands_per_task_raises(self, sample_config):
        with pytest.raises(ValueError, match="commands_per_task"):
            generate_search_array_script(
                sample_config, commands_per_task=0, array_max=259
            )

    def test_invalid_array_max_raises(self, sample_config):
        with pytest.raises(ValueError, match="array_max"):
            generate_search_array_script(
                sample_config, commands_per_task=50, array_max=-1
            )

    def test_missing_array_max_raises(self, sample_config):
        with pytest.raises(ValueError, match="array_max"):
            generate_search_array_script(
                sample_config, commands_per_task=50, array_max=None
            )


class TestResumeScript:
    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_resume" in generate_resume_script(sample_config)

    def test_uses_b_flag(self, sample_config):
        assert "-b " in generate_resume_script(sample_config)

    def test_uses_threads(self, sample_config):
        script = generate_resume_script(sample_config)
        assert "-t 16" in script
        assert "-a 8" in script

    def test_does_not_pass_search_program(self, sample_config):
        # The resume command line should not contain "-S diamond" — search
        # program is already configured in WorkingDirectory by prepare phase.
        assert "-S diamond" not in generate_resume_script(sample_config)

    def test_finds_working_directory(self, sample_config):
        assert "WorkingDirectory" in generate_resume_script(sample_config)


class TestCrossCutting:
    def _all_scripts(self, config):
        return [
            generate_prepare_script(config),
            _search_script(config),
            generate_resume_script(config),
        ]

    def test_all_scripts_have_set_euo_pipefail(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert "set -euo pipefail" in script

    def test_all_scripts_activate_conda(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert "conda activate convgeno" in script

    def test_raises_without_orthofinder_config(self):
        config = PipelineConfig(
            project_dir="/project", slurm=SlurmConfig(partition="hawkcpu")
        )
        # The orthofinder-missing check fires before array_max validation
        # in generate_search_array_script, so all three raise the same
        # "OrthoFinder" error even if we pass valid array_max.
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_prepare_script(config)
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_search_array_script(
                config, commands_per_task=50, array_max=99
            )
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_resume_script(config)
