"""Tests for convgeno.slurm.build_search_manifest (prepare-time glue)."""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.slurm.build_search_manifest import build_manifest, main


def _fake_workdir(tmp_path: Path, sizes: dict[int, int]) -> Path:
    work = tmp_path / "WorkingDirectory"
    work.mkdir()
    for sid, nbytes in sizes.items():
        # Species{id}.fa of the requested byte size (header + payload).
        (work / f"Species{sid}.fa").write_text(">g\n" + "M" * max(1, nbytes - 3))
    return work


def _search_commands(tmp_path: Path, pairs: list[tuple[int, int]], work: Path) -> Path:
    lines = [
        f"diamond blastp -d {work}/diamondDBSpecies{j} -q {work}/Species{i}.fa "
        f"-o {work}/Blast{i}_{j}.txt.gz -p 1"
        for i, j in pairs
    ]
    path = tmp_path / "run_search_commands.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


class TestBuildManifest:
    def test_end_to_end_partitions_all_commands(self, tmp_path):
        work = _fake_workdir(tmp_path, {0: 400, 1: 800})
        pairs = [(i, j) for i in range(2) for j in range(2)]  # 4 commands
        cmds = _search_commands(tmp_path, pairs, work)
        manifest_dir = tmp_path / "manifests"

        summary = build_manifest(cmds, work, manifest_dir, tasks=2)

        assert summary["commands"] == 4
        assert summary["tasks"] == 2
        # mem determinant is the largest Species*.fa (species 1, ~800 bytes).
        assert summary["mem_determinant_bytes"] == (work / "Species1.fa").stat().st_size
        files = sorted(manifest_dir.glob("search_task_*.txt"))
        assert [f.name for f in files] == ["search_task_0.txt", "search_task_1.txt"]
        seen = [ln for f in files for ln in f.read_text().splitlines()]
        assert len(seen) == 4
        assert all("Blast" in ln for ln in seen)

    def test_empty_commands_file_raises(self, tmp_path):
        work = _fake_workdir(tmp_path, {0: 100})
        empty = tmp_path / "empty.txt"
        empty.write_text("\n\n")
        with pytest.raises(ValueError, match="No search commands"):
            build_manifest(empty, work, tmp_path / "m", tasks=1)

    def test_main_success_returns_zero(self, tmp_path):
        work = _fake_workdir(tmp_path, {0: 100, 1: 100})
        cmds = _search_commands(tmp_path, [(0, 0), (0, 1), (1, 0), (1, 1)], work)
        manifest_dir = tmp_path / "m"
        rc = main(
            [
                "--search-commands", str(cmds),
                "--work-dir", str(work),
                "--manifest-dir", str(manifest_dir),
                "--tasks", "3",
            ]
        )
        assert rc == 0
        assert len(list(manifest_dir.glob("search_task_*.txt"))) == 3

    def test_main_missing_workdir_returns_one(self, tmp_path):
        cmds = _search_commands(tmp_path, [(0, 0)], tmp_path / "WorkingDirectory")
        rc = main(
            [
                "--search-commands", str(cmds),
                "--work-dir", str(tmp_path / "nonexistent"),
                "--manifest-dir", str(tmp_path / "m"),
                "--tasks", "1",
            ]
        )
        assert rc == 1
