"""Tests for convgeno.cli.clean_cmd (the `convgeno clean` step)."""

from __future__ import annotations

from pathlib import Path

from convgeno.cli.clean_cmd import run_clean
from convgeno.slurm.config import PipelineConfig, ProteomeInputConfig, SlurmConfig


def _write_raw(raw: Path) -> None:
    raw.mkdir(parents=True, exist_ok=True)
    # Unidentified simple headers: passed through verbatim, no duplicate-gene error.
    (raw / "Homo_sapiens.fa").write_text(">seq1\nMKV\n>seq2\nMKVAA\n", encoding="utf-8")
    (raw / "Felis_catus.fa").write_text(">seq1\nMK\n", encoding="utf-8")


class TestRunCleanExplicit:
    def test_explicit_dirs(self, tmp_path: Path):
        raw = tmp_path / "raw"
        _write_raw(raw)
        out = tmp_path / "clean"
        rc = run_clean(raw_dir=str(raw), out_dir=str(out))
        assert rc == 0
        assert (out / "Homo_sapiens.fa").exists()
        assert (out / "Felis_catus.fa").exists()

    def test_stats_json_written(self, tmp_path: Path):
        raw = tmp_path / "raw"
        _write_raw(raw)
        out = tmp_path / "clean"
        stats = tmp_path / "stats.json"
        rc = run_clean(raw_dir=str(raw), out_dir=str(out), stats_json=str(stats))
        assert rc == 0
        assert stats.exists()


class TestRunCleanFromConfig:
    def test_uses_config_proteome_input(self, tmp_path: Path):
        raw = tmp_path / "raw"
        _write_raw(raw)
        out = tmp_path / "clean"
        cfg = PipelineConfig(
            project_dir=str(tmp_path),
            slurm=SlurmConfig(partition="p"),
            proteome_input=ProteomeInputConfig(raw_dir=str(raw), cleaned_dir=str(out)),
        )
        cfg_path = tmp_path / "pipeline_config.yaml"
        cfg.save(cfg_path)
        rc = run_clean(config_path=str(cfg_path))
        assert rc == 0
        assert (out / "Homo_sapiens.fa").exists()

    def test_no_config_no_args_errors(self, tmp_path: Path, monkeypatch):
        monkeypatch.chdir(tmp_path)  # no pipeline config present here
        assert run_clean() == 1

    def test_config_without_raw_dir_errors(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir=str(tmp_path),
            slurm=SlurmConfig(partition="p"),
        )
        cfg_path = tmp_path / "pipeline_config.yaml"
        cfg.save(cfg_path)
        assert run_clean(config_path=str(cfg_path)) == 1
