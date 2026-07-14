"""
Tests for convgeno.cli.orthofinder_cmd multi-node functions

Run with:  pytest tests/test_cli_multinode.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from convgeno.cli.orthofinder_cmd import (
    run_generate_multinode,
    run_submit_multinode,
)
from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig
from convgeno.slurm.runtime import CondaRuntimeConfig


_TEST_RUNTIME = CondaRuntimeConfig(
    conda_module="miniforge3/24.3.0-0",
    conda_base=Path("/share/apps/miniforge3/24.3.0-0"),
    conda_env_prefix=Path("/home/user/.conda/envs/convgeno"),
)


@pytest.fixture()
def sample_config_file(tmp_path: Path) -> dict:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    for name in ["Sp1", "Sp2", "Sp3", "Sp4"]:
        (proteomes / f"{name}.fa").write_text(">gene1\nMKTLLIL\n")

    config = PipelineConfig(
        project_dir=str(tmp_path),
        conda_env="convgeno",
        slurm=SlurmConfig(partition="hawkcpu"),
        orthofinder=OrthoFinderConfig(
            input_dir=str(proteomes),
            output_dir=str(tmp_path / "results"),
        ),
        runtime=_TEST_RUNTIME,
    )
    config_path = tmp_path / "config.yaml"
    config.save(config_path)
    return {"config_path": str(config_path), "tmp_path": tmp_path}


def _mock_validation_pass(num_species: int = 1) -> MagicMock:
    validation = MagicMock()
    validation.is_valid.return_value = True
    validation.errors = []
    validation.warnings = []
    validation.fasta_files = [Path(f"species_{i}.faa") for i in range(num_species)]
    validation.summary.return_value = (
        f"Validation passed.\n  FASTA files found: {num_species}"
    )
    return validation


class TestRunGenerateMultinode:
    @patch(
        "convgeno.cli.orthofinder_cmd.validate_orthofinder_inputs",
        return_value=_mock_validation_pass(),
    )
    def test_creates_prepare_and_resume(self, _mock_validate, sample_config_file):
        # Search array is a v2 stub, so generate produces prepare + resume only.
        f = sample_config_file
        result = run_generate_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
        )
        assert set(result.keys()) == {"prepare", "resume"}
        for path in result.values():
            assert path.exists()
            assert path.read_text(encoding="utf-8").startswith("#!/bin/bash")


class TestRunSubmitMultinode:
    @patch(
        "convgeno.cli.orthofinder_cmd.validate_orthofinder_inputs",
        return_value=_mock_validation_pass(),
    )
    @patch("convgeno.cli.orthofinder_cmd.submit_sbatch")
    def test_refuses_chain_until_search_v2(
        self, mock_sbatch, _mock_validate, sample_config_file
    ):
        # The v2 search array is not implemented, so run_generate_multinode
        # omits "search" and the dependency chain must not be submitted.
        f = sample_config_file
        run_submit_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
            skip_confirm=True,
        )
        mock_sbatch.assert_not_called()
