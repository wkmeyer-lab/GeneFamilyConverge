"""
Tests for convgeno.slurm.config

Run with:  pytest tests/test_slurm_config.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.slurm.config import PipelineConfig, SlurmConfig, normalize_optional_account


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
        assert "#SBATCH --time=72:00:00" in lines
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
        assert loaded.slurm.time_limit == "72:00:00"
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


class TestOptionalAccountHandling:
    @pytest.mark.parametrize("account", [None, "", " ", "null", "None", "NULL"])
    def test_optional_account_values_are_omitted(self, account):
        cfg = SlurmConfig(partition="hawkcpu", account=account)
        lines = cfg.to_sbatch_lines()
        assert not any("--account" in line for line in lines)

    def test_valid_account_is_included(self):
        cfg = SlurmConfig(partition="hawkcpu", account="valid_alloc")
        lines = cfg.to_sbatch_lines()
        assert "#SBATCH --account=valid_alloc" in lines

    @pytest.mark.parametrize("account", [None, "", "null", "None"])
    def test_from_dict_normalizes_null_variants(self, account):
        data = {"partition": "hawkcpu"}
        if account is not None:
            data["account"] = account
        cfg = SlurmConfig.from_dict(data)
        assert cfg.account is None

    def test_from_dict_preserves_valid_account(self):
        cfg = SlurmConfig.from_dict({"partition": "hawkcpu", "account": "wym219"})
        assert cfg.account == "wym219"

    def test_to_dict_omits_normalized_null_account(self):
        cfg = SlurmConfig(partition="hawkcpu", account="null")
        lines = cfg.to_sbatch_lines()
        assert not any("--account" in line for line in lines)

    def test_save_load_roundtrip_null_account(self, tmp_path):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu", account=None),
        )
        path = tmp_path / "config.yaml"
        config.save(path)
        loaded = PipelineConfig.load(path)
        assert loaded.slurm.account is None
        assert not any("--account" in line for line in loaded.slurm.to_sbatch_lines())


class TestOpenFileLimitConfig:
    def test_default_is_none(self):
        cfg = SlurmConfig(partition="hawkcpu")
        assert cfg.open_file_limit is None

    def test_from_dict_with_open_file_limit(self):
        cfg = SlurmConfig.from_dict(
            {"partition": "hawkcpu", "open_file_limit": 8192}
        )
        assert cfg.open_file_limit == 8192

    def test_from_dict_without_open_file_limit(self):
        cfg = SlurmConfig.from_dict({"partition": "hawkcpu"})
        assert cfg.open_file_limit is None

    def test_roundtrip_with_open_file_limit(self, tmp_path):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu", open_file_limit=16384),
        )
        path = tmp_path / "config.yaml"
        config.save(path)
        loaded = PipelineConfig.load(path)
        assert loaded.slurm.open_file_limit == 16384

    def test_roundtrip_null_open_file_limit(self, tmp_path):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu", open_file_limit=None),
        )
        path = tmp_path / "config.yaml"
        config.save(path)
        loaded = PipelineConfig.load(path)
        assert loaded.slurm.open_file_limit is None
