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


class TestSearchArrayScript:
    def test_is_array_job(self, sample_config):
        assert "--array=" in generate_search_array_script(sample_config)

    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_search" in generate_search_array_script(sample_config)

    def test_reads_commands_file(self, sample_config):
        assert "diamond_commands.txt" in generate_search_array_script(sample_config)

    def test_uses_slurm_array_task_id(self, sample_config):
        assert "SLURM_ARRAY_TASK_ID" in generate_search_array_script(sample_config)

    def test_exits_if_no_work(self, sample_config):
        assert "exit 0" in generate_search_array_script(sample_config)

    def test_counts_failures(self, sample_config):
        assert "FAILED" in generate_search_array_script(sample_config)

    def test_uses_search_threads(self, sample_config):
        assert "--cpus-per-task=16" in generate_search_array_script(sample_config)

    def test_custom_commands_per_task(self, sample_config):
        script = generate_search_array_script(sample_config, commands_per_task=100)
        assert "COMMANDS_PER_TASK=100" in script


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
    def test_all_scripts_have_set_euo_pipefail(self, sample_config):
        for gen in (
            generate_prepare_script,
            generate_search_array_script,
            generate_resume_script,
        ):
            assert "set -euo pipefail" in gen(sample_config)

    def test_all_scripts_activate_conda(self, sample_config):
        for gen in (
            generate_prepare_script,
            generate_search_array_script,
            generate_resume_script,
        ):
            assert "conda activate convgeno" in gen(sample_config)

    def test_raises_without_orthofinder_config(self):
        config = PipelineConfig(
            project_dir="/project", slurm=SlurmConfig(partition="hawkcpu")
        )
        for gen in (
            generate_prepare_script,
            generate_search_array_script,
            generate_resume_script,
        ):
            with pytest.raises(ValueError, match="OrthoFinder"):
                gen(config)
