"""Tests for convgeno.cli.run_cmd (profile rendering + orchestrator script)."""

from __future__ import annotations

from pathlib import Path

from convgeno.cli.run_cmd import (
    _mem_to_mb,
    _orchestrator_script,
    _time_to_minutes,
    render_slurm_profile,
)
from convgeno.slurm.config import PipelineConfig, SlurmConfig
from convgeno.slurm.runtime import CondaRuntimeConfig


class TestMemParsing:
    def test_suffixes(self):
        assert _mem_to_mb("350400M") == 350400
        assert _mem_to_mb("120G") == 122880
        assert _mem_to_mb("2T") == 2097152
        assert _mem_to_mb("8000") == 8000  # bare number -> MB

    def test_none(self):
        assert _mem_to_mb(None) is None


class TestTimeParsing:
    def test_hms(self):
        assert _time_to_minutes("72:00:00") == 4320

    def test_days(self):
        assert _time_to_minutes("3-00:00:00") == 4320

    def test_minutes_seconds_round_up(self):
        assert _time_to_minutes("00:00:30") == 1  # sub-minute rounds up

    def test_none(self):
        assert _time_to_minutes(None) is None


class TestRenderProfile:
    def test_full(self):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(
                partition="hawkcpu",
                mem="120G",
                time_limit="72:00:00",
                account="myacct",
            ),
        )
        prof = render_slurm_profile(cfg, jobs=8)
        assert prof["executor"] == "slurm"
        assert prof["jobs"] == 8
        dr = prof["default-resources"]
        assert dr["slurm_partition"] == "hawkcpu"
        assert dr["slurm_account"] == "myacct"
        assert dr["mem_mb"] == 122880
        assert dr["runtime"] == 4320

    def test_omits_null_account(self):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu", account=None),
        )
        prof = render_slurm_profile(cfg, jobs=4)
        assert "slurm_account" not in prof.get("default-resources", {})


class TestOrchestratorScript:
    def _cfg(self, account):
        return PipelineConfig(
            project_dir="/proj",
            slurm=SlurmConfig(
                partition="hawkcpu", time_limit="72:00:00", account=account
            ),
            runtime=CondaRuntimeConfig(
                conda_module=None,
                conda_base=Path("/opt/conda"),
                conda_env_prefix=Path("/opt/conda/envs/convgeno"),
            ),
        )

    def test_has_account_and_activation(self):
        script = _orchestrator_script(
            self._cfg("myacct"),
            ["snakemake", "--snakefile", "/proj/workflow/Snakefile", "all"],
            Path("/proj"),
        )
        assert "#SBATCH --account=myacct" in script
        assert "#SBATCH --partition=hawkcpu" in script
        assert "conda activate" in script  # from the conda bootstrap
        assert "snakemake" in script
        assert 'cd "/proj"' in script

    def test_omits_account_when_none(self):
        script = _orchestrator_script(
            self._cfg(None),
            ["snakemake", "all"],
            Path("/proj"),
        )
        assert "--account" not in script
