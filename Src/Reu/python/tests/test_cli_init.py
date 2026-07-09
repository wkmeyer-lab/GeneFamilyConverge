"""
Tests for convgeno.cli.init_cmd

Run with:  pytest tests/test_cli_init.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from convgeno.cli.init_cmd import (
    _default_orthofinder_output_dir,
    _derive_orthofinder_threads,
    _parse_aligner_choice,
    run_init,
)
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


class TestDeriveOrthoFinderThreads:
    def test_threads_from_48_cpus(self):
        search_threads, analysis_threads = _derive_orthofinder_threads(48)

        assert search_threads == 48
        assert analysis_threads == 12

    def test_threads_from_hyperthreaded_detection(self):
        # Hyperthreading is handled upstream by detect_node_cpus:
        # 104 logical CPUs with ThreadsPerCore=2 gives 52 physical cores,
        # then reserves 4 cores, so init passes recommended_physical=48.
        search_threads, analysis_threads = _derive_orthofinder_threads(48)

        assert search_threads == 48
        assert analysis_threads == 12

    def test_threads_from_32_cpus(self):
        search_threads, analysis_threads = _derive_orthofinder_threads(32)

        assert search_threads == 32
        assert analysis_threads == 8

    def test_threads_from_8_cpus(self):
        search_threads, analysis_threads = _derive_orthofinder_threads(8)

        assert search_threads == 8
        assert analysis_threads == 2

    def test_threads_fallback_on_detection_failure(self):
        search_threads, analysis_threads = _derive_orthofinder_threads(0)

        assert search_threads == 16
        assert analysis_threads == 4


class TestParseAlignerChoice:
    def test_parse_aligner_default_empty(self):
        assert _parse_aligner_choice("") == "mafft"

    def test_parse_aligner_1(self):
        assert _parse_aligner_choice("1") == "mafft"

    def test_parse_aligner_2(self):
        assert _parse_aligner_choice("2") == "famsa"

    def test_parse_aligner_mafft_string(self):
        assert _parse_aligner_choice("mafft") == "mafft"

    def test_parse_aligner_famsa_string(self):
        assert _parse_aligner_choice("FAMSA") == "famsa"

    def test_parse_aligner_invalid_falls_back(self):
        assert _parse_aligner_choice("muscle") == "mafft"


class TestDefaultOrthoFinderOutputDir:
    def test_default_orthofinder_output_dir_uses_timestamp(self):
        path = _default_orthofinder_output_dir(
            Path("/project"),
            timestamp="20260629_181530",
        )

        assert path == Path(
            "/project/Data/processed/orthofinder_single_20260629_181530"
        )


class TestInitCreatesConfig:
    def test_happy_path_with_partitions(self, tmp_path: Path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        output = tmp_path / "pipeline_config.yaml"
        static_output_dir = tmp_path / "Data" / "processed" / "orthofinder"
        static_output_dir.mkdir(parents=True)
        (static_output_dir / ".gitkeep").write_text("", encoding="utf-8")

        inputs = iter([
            str(tmp_path),  # project directory
            "convgeno",     # conda env
            "1",            # select partition (hawkcpu)
            "",             # accept recommended cpus per task
            "",             # accept recommended memory
            "",             # accept default MSA aligner
            "72:00:00",     # time limit
            "",             # mail user (skip)
            "",             # account (skip)
            "",             # accept detected runtime defaults, if present
        ])

        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[HAWK, RAPIDS]),
            patch(
                "convgeno.cli.init_cmd.detect_node_cpus",
                return_value={
                    "min_cpus_per_node": 52,
                    "max_cpus_per_node": 52,
                    "recommended_cpus": 48,
                    "threads_per_core": 1,
                    "physical_cores": 52,
                    "recommended_physical": 48,
                    "node_count": 3,
                },
            ),
            patch(
                "convgeno.cli.init_cmd.detect_partition_memory",
                return_value={
                    "max_mem_per_cpu_mb": 7300,
                    "def_mem_per_cpu_mb": None,
                    "max_mem_per_node_mb": None,
                    "def_mem_per_node_mb": None,
                    "min_node_memory_mb": 380000,
                },
            ),
            patch(
                "convgeno.cli.init_cmd.detect_scratch_dir",
                return_value={
                    "scratch_base": None,
                    "is_ephemeral": False,
                    "method": "none",
                },
            ),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(output))

        assert output.exists()
        loaded = PipelineConfig.load(output)
        assert loaded.slurm.partition == "hawkcpu"
        assert loaded.slurm.cpus_per_task == 48
        assert loaded.slurm.mem == "350400M"
        assert loaded.slurm.mem_per_cpu is None
        assert loaded.conda_env == "convgeno"
        assert loaded.orthofinder is not None
        assert loaded.orthofinder.search_threads == 48
        assert loaded.orthofinder.analysis_threads == 12
        assert loaded.orthofinder.msa_program == "mafft"
        output_dir = Path(loaded.orthofinder.output_dir)
        assert output_dir.parent == tmp_path / "Data" / "processed"
        assert output_dir.name.startswith("orthofinder_single_")
        timestamp_suffix = output_dir.name.removeprefix("orthofinder_single_")
        assert len(timestamp_suffix) == 15
        assert timestamp_suffix[8] == "_"
        assert timestamp_suffix[:8].isdigit()
        assert timestamp_suffix[9:].isdigit()
        assert output_dir != static_output_dir
        assert not output_dir.exists()

    def test_no_partitions_detected(self, tmp_path: Path):
        output = tmp_path / "config.yaml"

        inputs = iter([
            "/some/project/path",  # project dir
            "myenv",               # conda env
            "gpu-partition",       # manually typed partition
            "32",                  # cpus
            "",                    # accept recommended memory
            "2",                   # famsa aligner
            "24:00:00",            # time limit
            "user@example.com",    # mail
            "myaccount",           # account
            "",                    # accept detected runtime defaults, if present
        ])

        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[]),
            patch(
                "convgeno.cli.init_cmd.detect_partition_memory",
                return_value={
                    "max_mem_per_cpu_mb": None,
                    "def_mem_per_cpu_mb": 5000,
                    "max_mem_per_node_mb": None,
                    "def_mem_per_node_mb": None,
                    "min_node_memory_mb": 256000,
                },
            ),
            patch(
                "convgeno.cli.init_cmd.detect_scratch_dir",
                return_value={
                    "scratch_base": None,
                    "is_ephemeral": False,
                    "method": "none",
                },
            ),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(output))

        loaded = PipelineConfig.load(output)
        assert loaded.slurm.partition == "gpu-partition"
        assert loaded.slurm.cpus_per_task == 32
        assert loaded.slurm.mem == "160000M"
        assert loaded.slurm.mail_user == "user@example.com"
        assert loaded.slurm.account == "myaccount"
        assert loaded.conda_env == "myenv"
        assert loaded.orthofinder is not None
        assert loaded.orthofinder.msa_program == "famsa"

    def test_memory_detection_failure_prompts_for_explicit_memory(
        self, tmp_path: Path, capsys
    ):
        output = tmp_path / "config.yaml"

        inputs = iter([
            "/some/project/path",  # project dir
            "convgeno",            # conda env
            "hawkcpu",             # manually typed partition
            "48",                  # cpus
            "350400M",             # explicit memory request
            "",                    # aligner
            "72:00:00",            # time limit
            "",                    # mail
            "",                    # account
            "",                    # accept detected runtime defaults, if present
        ])

        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[]),
            patch(
                "convgeno.cli.init_cmd.detect_partition_memory",
                return_value={
                    "max_mem_per_cpu_mb": None,
                    "def_mem_per_cpu_mb": None,
                    "max_mem_per_node_mb": None,
                    "def_mem_per_node_mb": None,
                    "min_node_memory_mb": None,
                },
            ),
            patch(
                "convgeno.cli.init_cmd.detect_scratch_dir",
                return_value={
                    "scratch_base": None,
                    "is_ephemeral": False,
                    "method": "none",
                },
            ),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(output))

        captured = capsys.readouterr()
        loaded = PipelineConfig.load(output)
        assert "Could not detect partition memory limits" in captured.out
        assert loaded.slurm.mem == "350400M"
        assert loaded.slurm.mem_per_cpu is None


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
            "120G",       # memory override
            "",           # aligner
            "12:00:00",   # time
            "",           # mail (skip)
            "",           # account (skip)
            "",           # accept detected runtime defaults, if present
        ])

        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[]),
            patch(
                "convgeno.cli.init_cmd.detect_partition_memory",
                return_value={
                    "max_mem_per_cpu_mb": None,
                    "def_mem_per_cpu_mb": None,
                    "max_mem_per_node_mb": None,
                    "def_mem_per_node_mb": None,
                    "min_node_memory_mb": 100000,
                },
            ),
            patch(
                "convgeno.cli.init_cmd.detect_scratch_dir",
                return_value={
                    "scratch_base": None,
                    "is_ephemeral": False,
                    "method": "none",
                },
            ),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(config_path))

        loaded = PipelineConfig.load(config_path)
        assert loaded.slurm.partition == "testpart"
        assert loaded.slurm.mem == "120G"


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
            "",             # memory
            "",             # aligner
            "72:00:00",     # time
            "",             # mail
            "",             # account
            "",             # accept detected runtime defaults, if present
        ])

        output = tmp_path / "config.yaml"
        with (
            patch("convgeno.cli.init_cmd.discover_partitions", return_value=[small]),
            patch(
                "convgeno.cli.init_cmd.detect_partition_memory",
                return_value={
                    "max_mem_per_cpu_mb": 7300,
                    "def_mem_per_cpu_mb": None,
                    "max_mem_per_node_mb": None,
                    "def_mem_per_node_mb": None,
                    "min_node_memory_mb": 16000,
                },
            ),
            patch(
                "convgeno.cli.init_cmd.detect_scratch_dir",
                return_value={
                    "scratch_base": None,
                    "is_ephemeral": False,
                    "method": "none",
                },
            ),
            patch("builtins.input", side_effect=inputs),
        ):
            run_init(output_path=str(output))

        captured = capsys.readouterr()
        assert "Warning" in captured.out
        assert output.exists()
