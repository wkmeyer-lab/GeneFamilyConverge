"""
Tests for convgeno.cli.orthofinder_cmd

Run with:  pytest tests/test_cli_orthofinder.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from convgeno.cli.orthofinder_cmd import (
    SubmitResult,
    is_invalid_account_error,
    run_generate,
    run_submit,
    strip_account_directive,
    submit_sbatch,
    validate_orthofinder_output_dir,
)
from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig, normalize_optional_account
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
        (proteomes / f"{name}.fa").write_text(f">gene1\nMKTLLIL\n")

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
    @patch("convgeno.cli.orthofinder_cmd._run_sbatch")
    def test_parses_job_id(self, mock_run, tmp_path):
        script = tmp_path / "job.sh"
        script.write_text("#!/bin/bash\n#SBATCH --partition=hawkcpu\necho hi\n")
        mock_run.return_value = MagicMock(
            returncode=0, stdout="Submitted batch job 98765432\n"
        )
        sr = submit_sbatch(script)
        assert sr.job_id == "98765432"
        assert sr.account_stripped is False

    @patch("convgeno.cli.orthofinder_cmd._run_sbatch")
    def test_raises_on_failure(self, mock_run, tmp_path):
        script = tmp_path / "job.sh"
        script.write_text("#!/bin/bash\necho hi\n")
        mock_run.return_value = MagicMock(
            returncode=1, stderr="sbatch: error: invalid partition"
        )
        with pytest.raises(RuntimeError, match="invalid partition"):
            submit_sbatch(script)

    @patch("convgeno.cli.orthofinder_cmd._run_sbatch")
    def test_raises_on_unparseable_output(self, mock_run, tmp_path):
        script = tmp_path / "job.sh"
        script.write_text("#!/bin/bash\necho hi\n")
        mock_run.return_value = MagicMock(
            returncode=0, stdout="some unexpected output\n"
        )
        with pytest.raises(RuntimeError, match="Could not parse"):
            submit_sbatch(script)


class TestRunSubmit:
    @patch(
        "convgeno.cli.orthofinder_cmd.submit_sbatch",
        return_value=SubmitResult(job_id="12345"),
    )
    @patch("builtins.input", return_value="y")
    def test_with_confirmation_yes(self, _mock_input, mock_sbatch, sample_config_file):
        f = sample_config_file
        run_submit(
            config_path=f["config_path"],
            script_path=str(f["tmp_path"] / "job.sh"),
            skip_confirm=False,
        )
        mock_sbatch.assert_called_once()

    @patch(
        "convgeno.cli.orthofinder_cmd.submit_sbatch",
        return_value=SubmitResult(job_id="12345"),
    )
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

    @patch(
        "convgeno.cli.orthofinder_cmd.submit_sbatch",
        return_value=SubmitResult(job_id="99999"),
    )
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


# ── Account normalization ─────────────────────────────────────────────

class TestNormalizeOptionalAccount:
    @pytest.mark.parametrize("raw", [None, "", " ", "null", "None", "NULL", " null "])
    def test_blank_and_null_variants_become_none(self, raw):
        assert normalize_optional_account(raw) is None

    def test_real_value_preserved(self):
        assert normalize_optional_account("wym219") == "wym219"

    def test_whitespace_stripped_from_real_value(self):
        assert normalize_optional_account("  my_alloc  ") == "my_alloc"


class TestAccountHelpers:
    def test_is_invalid_account_error_true(self):
        assert is_invalid_account_error(
            "Batch job submission failed: Invalid account or "
            "account/partition combination specified"
        )

    def test_is_invalid_account_error_false(self):
        assert not is_invalid_account_error("sbatch: error: invalid partition")

    def test_strip_account_directive_removes_line(self):
        script = (
            "#!/bin/bash\n"
            "#SBATCH --partition=hawkcpu\n"
            "#SBATCH --account=prm526\n"
            "#SBATCH --time=72:00:00\n"
            "echo hi\n"
        )
        stripped = strip_account_directive(script)
        assert "--account" not in stripped
        assert "#SBATCH --partition=hawkcpu" in stripped
        assert "#SBATCH --time=72:00:00" in stripped

    def test_strip_account_directive_noop_when_absent(self):
        script = "#!/bin/bash\n#SBATCH --partition=hawkcpu\necho hi\n"
        assert strip_account_directive(script) == script


# ── Account omission in generated scripts ─────────────────────────────

class TestAccountOmissionInSingleNodeScript:
    @pytest.mark.parametrize("account", [None, "", "null", "None"])
    def test_null_account_omits_sbatch_line(self, account, tmp_path):
        proteomes = tmp_path / "proteomes"
        proteomes.mkdir()
        (proteomes / "Sp1.fa").write_text(">g\nMK\n")
        config = PipelineConfig(
            project_dir=str(tmp_path),
            slurm=SlurmConfig(partition="hawkcpu", account=account),
            orthofinder=OrthoFinderConfig(
                input_dir=str(proteomes),
                output_dir=str(tmp_path / "results"),
            ),
            runtime=_TEST_RUNTIME,
        )
        from convgeno.slurm.script_generator import generate_orthofinder_script

        script = generate_orthofinder_script(config)
        assert "--account" not in script

    def test_valid_account_included(self, tmp_path):
        proteomes = tmp_path / "proteomes"
        proteomes.mkdir()
        (proteomes / "Sp1.fa").write_text(">g\nMK\n")
        config = PipelineConfig(
            project_dir=str(tmp_path),
            slurm=SlurmConfig(partition="hawkcpu", account="valid_alloc"),
            orthofinder=OrthoFinderConfig(
                input_dir=str(proteomes),
                output_dir=str(tmp_path / "results"),
            ),
            runtime=_TEST_RUNTIME,
        )
        from convgeno.slurm.script_generator import generate_orthofinder_script

        script = generate_orthofinder_script(config)
        assert "#SBATCH --account=valid_alloc" in script


# ── Submission fallback tests ─────────────────────────────────────────

class TestSubmitSbatchFallback:
    def _script_with_account(self, tmp_path):
        script = tmp_path / "job.sh"
        script.write_text(
            "#!/bin/bash\n"
            "#SBATCH --partition=hawkcpu\n"
            "#SBATCH --account=prm526\n"
            "echo hi\n"
        )
        return script

    def _script_without_account(self, tmp_path):
        script = tmp_path / "job.sh"
        script.write_text(
            "#!/bin/bash\n"
            "#SBATCH --partition=hawkcpu\n"
            "echo hi\n"
        )
        return script

    @patch("convgeno.cli.orthofinder_cmd._run_sbatch")
    def test_invalid_account_fallback_succeeds(self, mock_run, tmp_path):
        script = self._script_with_account(tmp_path)
        mock_run.side_effect = [
            MagicMock(
                returncode=1,
                stderr="Batch job submission failed: Invalid account or "
                       "account/partition combination specified",
            ),
            MagicMock(
                returncode=0,
                stdout="Submitted batch job 12345\n",
            ),
        ]
        sr = submit_sbatch(script)
        assert sr.job_id == "12345"
        assert sr.account_stripped is True
        assert mock_run.call_count == 2

    @patch("convgeno.cli.orthofinder_cmd._run_sbatch")
    def test_invalid_account_fallback_both_fail(self, mock_run, tmp_path):
        script = self._script_with_account(tmp_path)
        mock_run.side_effect = [
            MagicMock(
                returncode=1,
                stderr="Batch job submission failed: Invalid account",
            ),
            MagicMock(
                returncode=1,
                stderr="Batch job submission failed: some other error",
            ),
        ]
        with pytest.raises(RuntimeError, match="failed with and without"):
            submit_sbatch(script)

    @patch("convgeno.cli.orthofinder_cmd._run_sbatch")
    def test_no_account_cluster_requires_account(self, mock_run, tmp_path):
        script = self._script_without_account(tmp_path)
        mock_run.return_value = MagicMock(
            returncode=1,
            stderr="sbatch: error: You must specify a valid account",
        )
        with pytest.raises(RuntimeError, match="require a SLURM allocation account"):
            submit_sbatch(script)

    @patch("convgeno.cli.orthofinder_cmd._run_sbatch")
    def test_strip_account_flag_removes_line_upfront(self, mock_run, tmp_path):
        script = self._script_with_account(tmp_path)
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 67890\n",
        )
        sr = submit_sbatch(script, strip_account=True)
        assert sr.job_id == "67890"
        assert sr.account_stripped is True
        submitted_script = mock_run.call_args[0][0]
        content = Path(submitted_script).read_text()
        assert "--account" not in content
