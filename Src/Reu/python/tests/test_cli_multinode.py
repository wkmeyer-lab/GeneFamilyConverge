"""
Tests for convgeno.cli.orthofinder_cmd multi-node functions

Run with:  pytest tests/test_cli_multinode.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from convgeno.cli.orthofinder_cmd import (
    SubmitResult,
    run_generate_multinode,
    run_submit_multinode,
)
from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig


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
    def test_creates_three_scripts(self, _mock_validate, sample_config_file):
        f = sample_config_file
        result = run_generate_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
        )
        assert set(result.keys()) == {"prepare", "search", "resume"}
        for path in result.values():
            assert path.exists()
            assert path.read_text(encoding="utf-8").startswith("#!/bin/bash")


class TestRunSubmitMultinode:
    @patch(
        "convgeno.cli.orthofinder_cmd.validate_orthofinder_inputs",
        return_value=_mock_validation_pass(),
    )
    @patch(
        "convgeno.cli.orthofinder_cmd.submit_sbatch",
        side_effect=[
            SubmitResult(job_id="111"),
            SubmitResult(job_id="222"),
            SubmitResult(job_id="333"),
        ],
    )
    def test_chains_dependencies(
        self, mock_sbatch, _mock_validate, sample_config_file
    ):
        f = sample_config_file
        run_submit_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
            skip_confirm=True,
        )
        assert mock_sbatch.call_count == 3

        # First call: no dependency
        first_call = mock_sbatch.call_args_list[0]
        assert first_call.kwargs.get("dependency") is None

        # Second call: depends on Job A (111)
        second_call = mock_sbatch.call_args_list[1]
        assert second_call.kwargs.get("dependency") == "111"

        # Third call: depends on Job B (222)
        third_call = mock_sbatch.call_args_list[2]
        assert third_call.kwargs.get("dependency") == "222"

    @patch(
        "convgeno.cli.orthofinder_cmd.validate_orthofinder_inputs",
        return_value=_mock_validation_pass(),
    )
    @patch("convgeno.cli.orthofinder_cmd.submit_sbatch")
    @patch("builtins.input", return_value="n")
    def test_aborts_on_user_decline(
        self, _mock_input, mock_sbatch, _mock_validate, sample_config_file
    ):
        f = sample_config_file
        run_submit_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
            skip_confirm=False,
        )
        mock_sbatch.assert_not_called()


class TestArraySizeComputation:
    @patch(
        "convgeno.cli.orthofinder_cmd.validate_orthofinder_inputs",
        return_value=_mock_validation_pass(num_species=114),
    )
    def test_114_species_yields_array_0_to_259(
        self, _mock_validate, sample_config_file
    ):
        # 114 species → 114*114 = 12996 commands.
        # With commands_per_task=50, ceil(12996/50) = 260 tasks,
        # giving #SBATCH --array=0-259.
        f = sample_config_file
        result = run_generate_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
        )
        search_script = result["search"].read_text(encoding="utf-8")
        assert "#SBATCH --array=0-259" in search_script
        assert "#SBATCH --array=0-9999" not in search_script

    @patch(
        "convgeno.cli.orthofinder_cmd.validate_orthofinder_inputs",
        return_value=_mock_validation_pass(num_species=4),
    )
    def test_4_species_yields_small_array(
        self, _mock_validate, sample_config_file
    ):
        # 4 species → 16 commands → ceil(16/50) = 1 task → array=0-0.
        f = sample_config_file
        result = run_generate_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
        )
        search_script = result["search"].read_text(encoding="utf-8")
        assert "#SBATCH --array=0-0" in search_script
