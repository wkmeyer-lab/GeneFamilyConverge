"""
Tests for convgeno.slurm.discovery

Run with:  pytest tests/test_slurm_discovery.py -v
"""

from __future__ import annotations

from unittest.mock import Mock, patch

import pytest

from convgeno.slurm.discovery import (
    PartitionInfo,
    _parse_sinfo_line,
    detect_partition_memory,
    detect_node_cpus,
    detect_scratch_dir,
    discover_partitions,
    recommend_memory_mb,
)


class TestParseSinfoLine:
    def test_default_partition(self):
        info = _parse_sinfo_line("hawkcpu*         3-00:00:00     64     256000")
        assert info is not None
        assert info.name == "hawkcpu"
        assert info.is_default is True
        assert info.time_limit == "3-00:00:00"
        assert info.max_cpus_per_node == 64
        assert info.max_mem_mb_per_node == 256000

    def test_non_default(self):
        info = _parse_sinfo_line("rapids           3-00:00:00     128    512000")
        assert info is not None
        assert info.name == "rapids"
        assert info.is_default is False
        assert info.max_cpus_per_node == 128

    def test_memory_with_plus(self):
        info = _parse_sinfo_line("haswell          3-00:00:00     28     128000+")
        assert info is not None
        assert info.max_mem_mb_per_node == 128000

    def test_insufficient_fields(self):
        assert _parse_sinfo_line("hawkcpu 3-00:00:00") is None

    def test_empty_string(self):
        assert _parse_sinfo_line("") is None

    def test_nonnumeric_cpus(self):
        info = _parse_sinfo_line("badpart          3-00:00:00     N/A    256000")
        assert info is not None
        assert info.max_cpus_per_node == 0


MOCK_SINFO_OUTPUT = (
    "hawkcpu*         3-00:00:00     64     256000\n"
    "hawkcpu*         3-00:00:00     64     256000\n"
    "rapids           3-00:00:00     128    512000\n"
    "haswell          3-00:00:00     28     128000\n"
)


class TestDiscoverPartitions:
    @patch("convgeno.slurm.discovery.subprocess.run", side_effect=FileNotFoundError)
    def test_sinfo_not_available(self, mock_run):
        assert discover_partitions() == []

    @patch("convgeno.slurm.discovery.subprocess.run")
    def test_parses_multiline_output(self, mock_run):
        mock_run.return_value = Mock(returncode=0, stdout=MOCK_SINFO_OUTPUT)
        partitions = discover_partitions()
        assert len(partitions) == 3
        assert partitions[0].name == "hawkcpu"
        assert partitions[0].is_default is True
        remaining_names = {p.name for p in partitions[1:]}
        assert "rapids" in remaining_names
        assert "haswell" in remaining_names

    @patch("convgeno.slurm.discovery.subprocess.run")
    def test_nonzero_return_code(self, mock_run):
        mock_run.return_value = Mock(returncode=1, stdout="")
        assert discover_partitions() == []


class TestDetectNodeCpus:
    def test_typical_homogeneous_partition(self):
        with patch(
            "convgeno.slurm.discovery.subprocess.run",
            side_effect=[
                Mock(
                    returncode=0,
                    stdout="hawk-a001 52\nhawk-a002 52\nhawk-a003 52\n",
                ),
                Mock(returncode=0, stdout="NodeName=hawk-a001 ThreadsPerCore=1"),
            ],
        ):
            detected = detect_node_cpus("hawkcpu")

        assert detected["min_cpus_per_node"] == 52
        assert detected["max_cpus_per_node"] == 52
        assert detected["recommended_cpus"] == 48
        assert detected["threads_per_core"] == 1
        assert detected["physical_cores"] == 52
        assert detected["recommended_physical"] == 48
        assert detected["node_count"] == 3

    def test_heterogeneous_partition(self):
        with patch(
            "convgeno.slurm.discovery.subprocess.run",
            side_effect=[
                Mock(returncode=0, stdout="hawk-a001 52\nsol-f709 56\n"),
                Mock(returncode=0, stdout="NodeName=hawk-a001 ThreadsPerCore=1"),
            ],
        ):
            detected = detect_node_cpus("hawkcpu")

        assert detected["min_cpus_per_node"] == 52
        assert detected["max_cpus_per_node"] == 56
        assert detected["recommended_physical"] == 48

    def test_hyperthreaded_nodes(self):
        with patch(
            "convgeno.slurm.discovery.subprocess.run",
            side_effect=[
                Mock(returncode=0, stdout="hyper-a001 104\nhyper-a002 104\n"),
                Mock(returncode=0, stdout="NodeName=hyper-a001 ThreadsPerCore=2"),
            ],
        ):
            detected = detect_node_cpus("hyper")

        assert detected["physical_cores"] == 52
        assert detected["recommended_physical"] == 48

    def test_small_node(self):
        with patch(
            "convgeno.slurm.discovery.subprocess.run",
            side_effect=[
                Mock(returncode=0, stdout="small-a001 12\nsmall-a002 12\n"),
                Mock(returncode=0, stdout="NodeName=small-a001 ThreadsPerCore=1"),
            ],
        ):
            detected = detect_node_cpus("small")

        assert detected["recommended_physical"] == 8

    def test_very_small_node(self):
        with patch(
            "convgeno.slurm.discovery.subprocess.run",
            side_effect=[
                Mock(returncode=0, stdout="tiny-a001 8\n"),
                Mock(returncode=0, stdout="NodeName=tiny-a001 ThreadsPerCore=1"),
            ],
        ):
            detected = detect_node_cpus("tiny")

        assert detected["recommended_physical"] == 8

    @patch("convgeno.slurm.discovery.subprocess.run", side_effect=FileNotFoundError)
    def test_sinfo_fails(self, mock_run):
        detected = detect_node_cpus("missing")

        assert detected["min_cpus_per_node"] == 0
        assert detected["max_cpus_per_node"] == 0
        assert detected["recommended_cpus"] == 0
        assert detected["threads_per_core"] == 1
        assert detected["physical_cores"] == 0
        assert detected["recommended_physical"] == 16
        assert detected["node_count"] == 0

    def test_scontrol_fails_gracefully(self):
        with patch(
            "convgeno.slurm.discovery.subprocess.run",
            side_effect=[
                Mock(returncode=0, stdout="hawk-a001 52\nhawk-a002 52\n"),
                RuntimeError("scontrol unavailable"),
            ],
        ):
            detected = detect_node_cpus("hawkcpu")

        assert detected["min_cpus_per_node"] == 52
        assert detected["max_cpus_per_node"] == 52
        assert detected["threads_per_core"] == 1
        assert detected["physical_cores"] == 52
        assert detected["recommended_physical"] == 48


class TestDetectPartitionMemory:
    def test_detects_partition_and_min_node_memory(self):
        with patch(
            "convgeno.slurm.discovery.subprocess.run",
            side_effect=[
                Mock(
                    returncode=0,
                    stdout=(
                        "PartitionName=hawkcpu MaxMemPerCPU=7300 "
                        "DefMemPerCPU=UNLIMITED MaxMemPerNode=UNLIMITED "
                        "DefMemPerNode=0"
                    ),
                ),
                Mock(returncode=0, stdout="hawk-a001 380000\nhawk-a002 384000\n"),
            ],
        ):
            detected = detect_partition_memory("hawkcpu")

        assert detected["max_mem_per_cpu_mb"] == 7300
        assert detected["def_mem_per_cpu_mb"] is None
        assert detected["max_mem_per_node_mb"] is None
        assert detected["def_mem_per_node_mb"] is None
        assert detected["min_node_memory_mb"] == 380000

    @patch("convgeno.slurm.discovery.subprocess.run", side_effect=FileNotFoundError)
    def test_detect_partition_memory_fails_gracefully(self, mock_run):
        detected = detect_partition_memory("missing")

        assert detected["max_mem_per_cpu_mb"] is None
        assert detected["def_mem_per_cpu_mb"] is None
        assert detected["min_node_memory_mb"] is None


class TestRecommendMemoryMb:
    def test_uses_max_mem_per_cpu_first(self):
        assert recommend_memory_mb(48, 7300, 6000, 380000) == 350400

    def test_falls_back_to_def_mem_per_cpu(self):
        assert recommend_memory_mb(32, None, 5000, 380000) == 160000

    def test_falls_back_to_90_percent_min_node_memory(self):
        assert recommend_memory_mb(16, None, None, 380000) == 342000

    def test_raises_when_nothing_available(self):
        with pytest.raises(ValueError, match="Could not detect partition memory"):
            recommend_memory_mb(16, None, None, None)


class TestDetectScratchDir:
    def test_scratch_env_var(self, tmp_path):
        with patch.dict(
            "os.environ",
            {"SCRATCH": str(tmp_path), "USER": "testuser"},
            clear=True,
        ):
            detected = detect_scratch_dir()

        assert detected["scratch_base"] == str(tmp_path)
        assert detected["method"] == "env_var"

    def test_ceph_scratch_found(self):
        username = "testuser"
        scratch_path = f"/share/ceph/scratch/{username}"

        with (
            patch.dict("os.environ", {"USER": username}, clear=True),
            patch(
                "convgeno.slurm.discovery.os.path.isdir",
                side_effect=lambda path: path == scratch_path,
            ),
            patch("convgeno.slurm.discovery.os.access", return_value=True),
        ):
            detected = detect_scratch_dir()

        assert detected["scratch_base"] == scratch_path
        assert detected["is_ephemeral"] is False
        assert detected["method"] == "path_probe"

    def test_local_scratch_is_ephemeral(self):
        with (
            patch.dict("os.environ", {"USER": "testuser"}, clear=True),
            patch(
                "convgeno.slurm.discovery.os.path.isdir",
                side_effect=lambda path: path == "/local/scratch",
            ),
            patch("convgeno.slurm.discovery.os.access", return_value=True),
        ):
            detected = detect_scratch_dir()

        assert detected["scratch_base"] == "/local/scratch"
        assert detected["is_ephemeral"] is True
        assert detected["method"] == "path_probe"

    def test_no_scratch_found(self):
        with (
            patch.dict("os.environ", {"USER": "testuser"}, clear=True),
            patch("convgeno.slurm.discovery.os.path.isdir", return_value=False),
            patch("convgeno.slurm.discovery.os.access", return_value=False),
        ):
            detected = detect_scratch_dir()

        assert detected["scratch_base"] is None
        assert detected["method"] == "none"

    def test_scratch_not_writable_skipped(self):
        with (
            patch.dict("os.environ", {"USER": "testuser"}, clear=True),
            patch("convgeno.slurm.discovery.os.path.isdir", return_value=True),
            patch("convgeno.slurm.discovery.os.access", return_value=False),
        ):
            detected = detect_scratch_dir()

        assert detected["scratch_base"] is None

    def test_filesystem_error_handled(self):
        with (
            patch.dict("os.environ", {"USER": "testuser"}, clear=True),
            patch(
                "convgeno.slurm.discovery.os.path.isdir",
                side_effect=OSError("filesystem unavailable"),
            ),
        ):
            detected = detect_scratch_dir()

        assert detected["scratch_base"] is None
        assert detected["method"] == "none"


class TestPartitionInfoStr:
    def test_default(self):
        info = PartitionInfo(
            name="hawkcpu",
            is_default=True,
            time_limit="3-00:00:00",
            max_cpus_per_node=64,
            max_mem_mb_per_node=256000,
        )
        s = str(info)
        assert "hawkcpu" in s
        assert "default" in s
        assert "3-00:00:00" in s
        assert "64" in s

    def test_non_default(self):
        info = PartitionInfo(
            name="rapids",
            is_default=False,
            time_limit="3-00:00:00",
            max_cpus_per_node=128,
            max_mem_mb_per_node=512000,
        )
        s = str(info)
        assert "rapids" in s
        assert "default" not in s
