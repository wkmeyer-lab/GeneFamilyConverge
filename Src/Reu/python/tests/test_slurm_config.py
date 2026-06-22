"""
Tests for convgeno.slurm.config

Run with:  pytest tests/test_slurm_config.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.slurm.config import PipelineConfig, SlurmConfig


class TestSlurmConfigFromDict:
    def test_valid_dict(self):
        cfg = SlurmConfig.from_dict(
            {"partition": "hawkcpu", "cpus_per_task": 32, "time_limit": "24:00:00"}
        )
        assert cfg.partition == "hawkcpu"
        assert cfg.cpus_per_task == 32
        assert cfg.time_limit == "24:00:00"
        assert cfg.nodes == 1
        assert cfg.ntasks == 1
        assert cfg.mem_per_cpu == "4G"
        assert cfg.account is None
        assert cfg.mail_user is None
        assert cfg.mail_type == "END,FAIL"
        assert cfg.extra_sbatch_args == []

    def test_missing_partition_raises(self):
        with pytest.raises(ValueError, match="partition"):
            SlurmConfig.from_dict({"cpus_per_task": 16})

    def test_ignores_unknown_keys(self):
        cfg = SlurmConfig.from_dict(
            {"partition": "rapids", "some_future_field": "value123"}
        )
        assert cfg.partition == "rapids"


class TestToSbatchLines:
    def test_basic(self):
        cfg = SlurmConfig(partition="hawkcpu")
        lines = cfg.to_sbatch_lines()
        assert isinstance(lines, list)
        assert "#SBATCH --partition=hawkcpu" in lines
        assert "#SBATCH --time=48:00:00" in lines
        assert "#SBATCH --cpus-per-task=16" in lines
        assert "#SBATCH --nodes=1" in lines
        assert not any("--account" in line for line in lines)
        assert not any("--mail-user" in line for line in lines)

    def test_with_extra_args(self):
        cfg = SlurmConfig(
            partition="rapids",
            extra_sbatch_args=["--qos=high", "--reservation=myres"],
        )
        lines = cfg.to_sbatch_lines()
        assert "#SBATCH --qos=high" in lines
        assert "#SBATCH --reservation=myres" in lines

    def test_with_optional_fields(self):
        cfg = SlurmConfig(
            partition="hawkcpu",
            account="wym219",
            mail_user="abc@lehigh.edu",
        )
        lines = cfg.to_sbatch_lines()
        assert "#SBATCH --account=wym219" in lines
        assert "#SBATCH --mail-user=abc@lehigh.edu" in lines


class TestPipelineConfigRoundtrip:
    def test_save_and_load(self, tmp_path: Path):
        original = PipelineConfig(
            project_dir="/share/ceph/wym219group/shared/projects/GeneFamilyConverge",
            conda_env="convgeno",
            slurm=SlurmConfig(partition="hawkcpu", cpus_per_task=32),
        )
        config_path = tmp_path / "config.yaml"
        original.save(config_path)
        loaded = PipelineConfig.load(config_path)

        assert loaded.project_dir == original.project_dir
        assert loaded.conda_env == original.conda_env
        assert loaded.slurm.partition == "hawkcpu"
        assert loaded.slurm.cpus_per_task == 32
        assert loaded.slurm.time_limit == "48:00:00"
        assert loaded.slurm.nodes == 1

    def test_load_missing_file(self):
        with pytest.raises(FileNotFoundError):
            PipelineConfig.load("/nonexistent/path/config.yaml")


class TestToDictExcludesNone:
    def test_none_fields_omitted(self):
        cfg = SlurmConfig(partition="hawkcpu")
        d = cfg.to_dict()
        assert "account" not in d
        assert "mail_user" not in d
        assert "partition" in d
