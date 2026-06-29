"""
Tests for convgeno.cli.init_cmd

Run with:  pytest tests/test_cli_init.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from convgeno.cli.init_cmd import run_init
from convgeno.slurm.config import PipelineConfig
from convgeno.slurm.discovery import PartitionInfo


HAWK = PartitionInfo(
    name="hawkcpu",
    is_default=True,
    time_limit="3-00:00:00",
    max_cpus_per_node=64,
    max_mem_mb_per_node=256000,
)
RAPIDS = PartitionInfo(
    name="rapids",
    is_default=False,
    time_limit="3-00:00:00",
    max_cpus_per_node=128,
    max_mem_mb_per_node=512000,
)


class TestInitCreatesConfig:
    def test_happy_path_with_partitions(self, tmp_path: Path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        output = tmp_path / "pipeline_config.yaml"

        inputs = iter([
            str(tmp_path),  # project directory
            "convgeno",     # conda env
            "1",            # select partition (hawkcpu)
            "16",           # cpus per task
            "72:00:00",     # time limit
            "4G",           # memory per cpu
            "",             # mail user (skip)
            "",             # account (skip)
            "",             # open_file_limit (default 8192)
            "",             # accept detected runtime defaults, if present
        ])

        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[HAWK, RAPIDS]),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(output))

        assert output.exists()
        loaded = PipelineConfig.load(output)
        assert loaded.slurm.partition == "hawkcpu"
        assert loaded.slurm.cpus_per_task == 16
        assert loaded.conda_env == "convgeno"

    def test_no_partitions_detected(self, tmp_path: Path):
        output = tmp_path / "config.yaml"

        inputs = iter([
            "/some/project/path",  # project dir
            "myenv",               # conda env
            "gpu-partition",       # manually typed partition
            "32",                  # cpus
            "24:00:00",            # time limit
            "8G",                  # memory
            "user@example.com",    # mail
            "myaccount",           # account
            "",                    # open_file_limit (default 8192)
            "",                    # accept detected runtime defaults, if present
        ])

        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[]),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(output))

        loaded = PipelineConfig.load(output)
        assert loaded.slurm.partition == "gpu-partition"
        assert loaded.slurm.cpus_per_task == 32
        assert loaded.slurm.mail_user == "user@example.com"
        assert loaded.slurm.account == "myaccount"
        assert loaded.conda_env == "myenv"


class TestInitOverwriteBehaviour:
    def test_aborts_if_user_declines(self, tmp_path: Path):
        config_path = tmp_path / "config.yaml"
        config_path.write_text("existing content", encoding="utf-8")

        with patch("builtins.input", return_value="n"):
            run_init(output_path=str(config_path))

        assert config_path.read_text(encoding="utf-8") == "existing content"

    def test_overwrites_if_user_confirms(self, tmp_path: Path):
        config_path = tmp_path / "config.yaml"
        config_path.write_text("old", encoding="utf-8")

        inputs = iter([
            "y",          # confirm overwrite
            "/project",   # project dir
            "convgeno",   # conda env
            "testpart",   # partition (no sinfo)
            "8",          # cpus
            "12:00:00",   # time
            "2G",         # mem
            "",           # mail (skip)
            "",           # account (skip)
            "",           # open_file_limit (default 8192)
            "",           # accept detected runtime defaults, if present
        ])

        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[]),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(config_path))

        loaded = PipelineConfig.load(config_path)
        assert loaded.slurm.partition == "testpart"


class TestInitWarnings:
    def test_warns_cpu_exceeds_partition(self, tmp_path: Path, capsys):
        small = PartitionInfo(
            name="small",
            is_default=True,
            time_limit="1:00:00",
            max_cpus_per_node=8,
            max_mem_mb_per_node=16000,
        )

        inputs = iter([
            str(tmp_path),  # project dir
            "convgeno",     # env
            "1",            # select partition
            "32",           # cpus (exceeds 8)
            "72:00:00",     # time
            "4G",           # mem
            "",             # mail
            "",             # account
            "",             # open_file_limit (default 8192)
            "",             # accept detected runtime defaults, if present
        ])

        output = tmp_path / "config.yaml"
        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[small]),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(output))

        captured = capsys.readouterr()
        assert "Warning" in captured.out
        assert output.exists()
