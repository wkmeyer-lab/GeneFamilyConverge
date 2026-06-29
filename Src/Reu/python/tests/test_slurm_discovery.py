"""
Tests for convgeno.slurm.discovery

Run with:  pytest tests/test_slurm_discovery.py -v
"""

from __future__ import annotations

from subprocess import CompletedProcess
from unittest.mock import patch

import pytest

from convgeno.slurm.discovery import (
    PartitionInfo,
    _parse_sinfo_line,
    discover_partitions,
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
        mock_run.return_value = CompletedProcess(
            args=[], returncode=0, stdout=MOCK_SINFO_OUTPUT,
        )
        partitions = discover_partitions()
        assert len(partitions) == 3
        assert partitions[0].name == "hawkcpu"
        assert partitions[0].is_default is True
        remaining_names = {p.name for p in partitions[1:]}
        assert "rapids" in remaining_names
        assert "haswell" in remaining_names

    @patch("convgeno.slurm.discovery.subprocess.run")
    def test_nonzero_return_code(self, mock_run):
        mock_run.return_value = CompletedProcess(
            args=[], returncode=1, stdout="",
        )
        assert discover_partitions() == []


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
