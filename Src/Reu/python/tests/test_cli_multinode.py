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
    def test_creates_all_three_scripts(self, _mock_validate, sample_config_file):
        # v2 generates prepare + search array + resume.
        f = sample_config_file
        result = run_generate_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
        )
        assert set(result.keys()) == {"prepare", "search", "resume"}
        for path in result.values():
            assert path.exists()
            assert path.read_text(encoding="utf-8").startswith("#!/bin/bash")
        # The search script is a SLURM array; the prepare builds the manifest.
        assert "--array=" in result["search"].read_text(encoding="utf-8")
        assert (
            "build_search_manifest"
            in result["prepare"].read_text(encoding="utf-8")
        )


class TestRunSubmitMultinode:
    @patch(
        "convgeno.cli.orthofinder_cmd.validate_orthofinder_inputs",
        return_value=_mock_validation_pass(),
    )
    @patch(
        "convgeno.cli.orthofinder_cmd.submit_sbatch",
        return_value=SubmitResult(job_id="123", account_stripped=False),
    )
    def test_submits_three_job_chain(
        self, mock_sbatch, _mock_validate, sample_config_file
    ):
        # v2 submits prepare -> search -> resume as a dependency chain.
        f = sample_config_file
        run_submit_multinode(
            config_path=f["config_path"],
            script_dir=str(f["tmp_path"] / "scripts"),
            skip_confirm=True,
        )
        assert mock_sbatch.call_count == 3
        # search + resume are submitted with a dependency on the prior job.
        search_call, resume_call = mock_sbatch.call_args_list[1:3]
        assert search_call.kwargs.get("dependency") == "123"
        assert resume_call.kwargs.get("dependency") == "123"
