"""
Tests for convgeno.cli.orthofinder_cmd

Run with:  pytest tests/test_cli_orthofinder.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from convgeno.cli.orthofinder_cmd import (
    run_generate,
    run_submit,
    submit_sbatch,
    validate_orthofinder_output_dir,
)
from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig


@pytest.fixture()
def sample_config_file(tmp_path: Path) -> dict:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    for name in ["Sp1", "Sp2", "Sp3", "Sp4"]:
        (proteomes / f"{name}.fa").write_text(f">gene1\nMKTLLIL\n")

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
    return {
        "config_path": str(config_path),
        "tmp_path": tmp_path,
        "config": config,
    }


class TestRunGenerate:
    def test_creates_script(self, sample_config_file):
        f = sample_config_file
        path = run_generate(
            config_path=f["config_path"],
            script_path=str(f["tmp_path"] / "job.sh"),
        )
        assert path.exists()
        content = path.read_text(encoding="utf-8")
        assert "#!/bin/bash" in content
        assert "orthofinder" in content
        assert "#SBATCH --partition=hawkcpu" in content

    def test_raises_if_config_missing(self):
        with pytest.raises(FileNotFoundError):
            run_generate(
                config_path="/nonexistent/config.yaml",
                script_path="/tmp/job.sh",
            )


class TestSubmitSbatch:
    @patch("convgeno.cli.orthofinder_cmd.subprocess.run")
    def test_parses_job_id(self, mock_run):
        mock_run.return_value = MagicMock(
            returncode=0, stdout="Submitted batch job 98765432\n"
        )
        assert submit_sbatch(Path("/fake/script.sh")) == "98765432"

    @patch("convgeno.cli.orthofinder_cmd.subprocess.run")
    def test_raises_on_failure(self, mock_run):
        mock_run.return_value = MagicMock(
            returncode=1, stderr="sbatch: error: invalid partition"
        )
        with pytest.raises(RuntimeError, match="invalid partition"):
            submit_sbatch(Path("/fake/script.sh"))

    @patch("convgeno.cli.orthofinder_cmd.subprocess.run")
    def test_raises_on_unparseable_output(self, mock_run):
        mock_run.return_value = MagicMock(
            returncode=0, stdout="some unexpected output\n"
        )
        with pytest.raises(RuntimeError, match="Could not parse"):
            submit_sbatch(Path("/fake/script.sh"))


class TestRunSubmit:
    @patch("convgeno.cli.orthofinder_cmd.submit_sbatch", return_value="12345")
    @patch("builtins.input", return_value="y")
    def test_with_confirmation_yes(self, _mock_input, mock_sbatch, sample_config_file):
        f = sample_config_file
        run_submit(
            config_path=f["config_path"],
            script_path=str(f["tmp_path"] / "job.sh"),
            skip_confirm=False,
        )
        mock_sbatch.assert_called_once()

    @patch("convgeno.cli.orthofinder_cmd.submit_sbatch", return_value="12345")
    @patch("builtins.input", return_value="n")
    def test_with_confirmation_no(self, _mock_input, mock_sbatch, sample_config_file):
        f = sample_config_file
        script = f["tmp_path"] / "job.sh"
        run_submit(
            config_path=f["config_path"],
            script_path=str(script),
            skip_confirm=False,
        )
        mock_sbatch.assert_not_called()
        assert script.exists()

    @patch("convgeno.cli.orthofinder_cmd.submit_sbatch", return_value="99999")
    def test_skip_confirm(self, mock_sbatch, sample_config_file):
        f = sample_config_file
        run_submit(
            config_path=f["config_path"],
            script_path=str(f["tmp_path"] / "job.sh"),
            skip_confirm=True,
        )
        mock_sbatch.assert_called_once()

    @patch(
        "convgeno.cli.orthofinder_cmd.submit_sbatch",
        side_effect=FileNotFoundError,
    )
    def test_sbatch_not_found(self, _mock_sbatch, sample_config_file):
        f = sample_config_file
        with pytest.raises(SystemExit):
            run_submit(
                config_path=f["config_path"],
                script_path=str(f["tmp_path"] / "job.sh"),
                skip_confirm=True,
            )


class TestValidateOrthofinderOutputDir:
    def test_rejects_existing_dir(self, tmp_path: Path):
        existing = tmp_path / "already_there"
        existing.mkdir()
        with pytest.raises(ValueError, match="already exists"):
            validate_orthofinder_output_dir(str(existing))

    def test_rejects_existing_file(self, tmp_path: Path):
        existing = tmp_path / "a_file"
        existing.write_text("x")
        with pytest.raises(ValueError, match="already exists"):
            validate_orthofinder_output_dir(str(existing))

    def test_rejects_single_quote_characters(self, tmp_path: Path):
        with pytest.raises(ValueError, match="quote characters"):
            validate_orthofinder_output_dir(f"'{tmp_path}/out'")

    def test_rejects_double_quote_characters(self, tmp_path: Path):
        with pytest.raises(ValueError, match="quote characters"):
            validate_orthofinder_output_dir(f'"{tmp_path}/out"')

    def test_accepts_nonexistent_dir(self, tmp_path: Path):
        fresh = tmp_path / "fresh_run"
        validate_orthofinder_output_dir(str(fresh))
        assert not fresh.exists()

    def test_creates_parent_directory(self, tmp_path: Path):
        nested = tmp_path / "parent" / "child" / "run"
        validate_orthofinder_output_dir(str(nested))
        assert nested.parent.is_dir()
        assert not nested.exists()


class TestRunGenerateOutputDirValidation:
    def test_fails_if_output_dir_already_exists(self, tmp_path: Path):
        # Build a config where output_dir already exists on disk.
        proteomes = tmp_path / "proteomes"
        proteomes.mkdir()
        for name in ["Sp1", "Sp2", "Sp3", "Sp4"]:
            (proteomes / f"{name}.fa").write_text(">g\nMK\n")
        existing_output = tmp_path / "already_exists"
        existing_output.mkdir()

        config = PipelineConfig(
            project_dir=str(tmp_path),
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(
                input_dir=str(proteomes),
                output_dir=str(existing_output),
            ),
        )
        config_path = tmp_path / "config.yaml"
        config.save(config_path)

        with pytest.raises(SystemExit):
            run_generate(
                config_path=str(config_path),
                script_path=str(tmp_path / "job.sh"),
            )
