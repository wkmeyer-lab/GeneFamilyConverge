"""
Tests for convgeno.external.config.OrthoFinderConfig

Run with:  pytest tests/test_orthofinder_config.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig


class TestOrthoFinderConfigFromDict:
    def test_valid_dict(self):
        cfg = OrthoFinderConfig.from_dict(
            {"input_dir": "/data/proteomes", "output_dir": "/data/results",
             "search_threads": 32}
        )
        assert cfg.input_dir == "/data/proteomes"
        assert cfg.output_dir == "/data/results"
        assert cfg.search_threads == 32
        assert cfg.analysis_threads == 8
        assert cfg.sequence_search == "diamond"
        assert cfg.msa_program == "mafft"
        assert cfg.tree_program == "fasttree"
        assert cfg.extra_args == []

    def test_missing_input_dir(self):
        with pytest.raises(ValueError, match="input_dir"):
            OrthoFinderConfig.from_dict({"output_dir": "/data/results"})

    def test_missing_output_dir(self):
        with pytest.raises(ValueError, match="output_dir"):
            OrthoFinderConfig.from_dict({"input_dir": "/data/proteomes"})

    def test_ignores_unknown_keys(self):
        cfg = OrthoFinderConfig.from_dict(
            {"input_dir": "/in", "output_dir": "/out", "future_setting": "xyz"}
        )
        assert cfg.input_dir == "/in"


class TestToCommandArgs:
    def test_default_args(self):
        cfg = OrthoFinderConfig(input_dir="/in", output_dir="/out")
        args = cfg.to_command_args()
        assert isinstance(args, list)
        assert args[args.index("-f") + 1] == "/in"
        assert args[args.index("-o") + 1] == "/out"
        assert args[args.index("-t") + 1] == "16"
        assert args[args.index("-a") + 1] == "8"
        assert args[args.index("-S") + 1] == "diamond"

    def test_with_extra_args(self):
        cfg = OrthoFinderConfig(
            input_dir="/in", output_dir="/out", extra_args=["--fewer-files"]
        )
        assert "--fewer-files" in cfg.to_command_args()


class TestToDictRoundtrip:
    def test_roundtrip(self):
        original = OrthoFinderConfig(
            input_dir="/in", output_dir="/out", search_threads=24
        )
        restored = OrthoFinderConfig.from_dict(original.to_dict())
        assert restored.input_dir == original.input_dir
        assert restored.output_dir == original.output_dir
        assert restored.search_threads == original.search_threads
        assert restored.analysis_threads == original.analysis_threads
        assert restored.sequence_search == original.sequence_search
        assert restored.msa_program == original.msa_program
        assert restored.tree_program == original.tree_program
        assert restored.extra_args == original.extra_args


class TestPipelineConfigWithOrthoFinder:
    def test_roundtrip_with_orthofinder(self, tmp_path: Path):
        of_config = OrthoFinderConfig(
            input_dir="/data/proteomes", output_dir="/data/results"
        )
        slurm = SlurmConfig(partition="hawkcpu")
        config = PipelineConfig(
            project_dir="/project", slurm=slurm, orthofinder=of_config
        )
        config_path = tmp_path / "config.yaml"
        config.save(config_path)
        loaded = PipelineConfig.load(config_path)

        assert loaded.orthofinder is not None
        assert loaded.orthofinder.input_dir == "/data/proteomes"
        assert loaded.orthofinder.output_dir == "/data/results"
        assert loaded.orthofinder.sequence_search == "diamond"

    def test_without_orthofinder_still_works(self, tmp_path: Path):
        config = PipelineConfig(
            project_dir="/project", slurm=SlurmConfig(partition="hawkcpu")
        )
        config_path = tmp_path / "config.yaml"
        config.save(config_path)
        loaded = PipelineConfig.load(config_path)

        assert loaded.orthofinder is None
        assert loaded.slurm.partition == "hawkcpu"
