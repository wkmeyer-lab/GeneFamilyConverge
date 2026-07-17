"""Tests for convgeno.slurm.verify_search_complete (resume completeness gate)."""

from __future__ import annotations

from pathlib import Path

from convgeno.slurm.verify_search_complete import (
    find_missing_blast_files,
    main,
    map_pairs_to_task_ids,
    verify,
)


def _workdir(tmp_path: Path, ids: list[int], present: list[tuple[int, int]]) -> Path:
    """Build a WorkingDirectory with SpeciesIDs.txt and the given Blast files."""
    work = tmp_path / "WorkingDirectory"
    work.mkdir()
    (work / "SpeciesIDs.txt").write_text(
        "\n".join(f"{i}: species{i}.fa" for i in ids) + "\n"
    )
    for i, j in present:
        (work / f"Blast{i}_{j}.txt.gz").write_text("data")
    return work


def _manifests(tmp_path: Path, buckets: dict[int, list[tuple[int, int]]]) -> Path:
    """Build a manifest dir; buckets maps task_id -> list of (i, j) pairs."""
    mdir = tmp_path / "manifests"
    mdir.mkdir()
    for task_id, pairs in buckets.items():
        lines = [
            f"diamond blastp -d /w/diamondDBSpecies{j} -q /w/Species{i}.fa "
            f"-o /w/Blast{i}_{j}.txt.gz -p 1"
            for i, j in pairs
        ]
        (mdir / f"search_task_{task_id}.txt").write_text("\n".join(lines) + "\n")
    return mdir


class TestFindMissing:
    def test_all_present(self, tmp_path):
        pairs = [(i, j) for i in range(2) for j in range(2)]
        work = _workdir(tmp_path, [0, 1], pairs)
        assert find_missing_blast_files(work, [0, 1]) == []

    def test_reports_missing_pair(self, tmp_path):
        work = _workdir(tmp_path, [0, 1], [(0, 0), (0, 1), (1, 0)])  # no (1,1)
        assert find_missing_blast_files(work, [0, 1]) == [(1, 1)]

    def test_empty_file_counts_as_missing(self, tmp_path):
        work = _workdir(tmp_path, [0], [])
        (work / "Blast0_0.txt.gz").write_text("")  # zero bytes
        assert find_missing_blast_files(work, [0]) == [(0, 0)]


class TestMapPairsToTaskIds:
    def test_maps_missing_pair_to_owning_task(self, tmp_path):
        mdir = _manifests(tmp_path, {0: [(0, 0), (1, 1)], 1: [(0, 1), (1, 0)]})
        assert map_pairs_to_task_ids(mdir, [(1, 0)]) == [1]
        assert map_pairs_to_task_ids(mdir, [(0, 0), (0, 1)]) == [0, 1]

    def test_no_substring_collision(self, tmp_path):
        # Blast0_1 must NOT match the task holding Blast0_10.
        mdir = _manifests(tmp_path, {0: [(0, 1)], 1: [(0, 10)]})
        assert map_pairs_to_task_ids(mdir, [(0, 1)]) == [0]
        assert map_pairs_to_task_ids(mdir, [(0, 10)]) == [1]

    def test_missing_manifest_dir_returns_empty(self, tmp_path):
        assert map_pairs_to_task_ids(tmp_path / "nope", [(0, 0)]) == []


class TestVerify:
    def test_ok_when_complete(self, tmp_path):
        pairs = [(i, j) for i in range(2) for j in range(2)]
        work = _workdir(tmp_path, [0, 1], pairs)
        mdir = _manifests(tmp_path, {0: pairs})
        ok, report = verify(work, mdir, "orthofinder_search.sh")
        assert ok is True
        assert "4/4" in report

    def test_reports_and_fails_when_incomplete(self, tmp_path):
        work = _workdir(tmp_path, [0, 1], [(0, 0), (0, 1), (1, 0)])
        mdir = _manifests(tmp_path, {0: [(0, 0), (0, 1)], 1: [(1, 0), (1, 1)]})
        ok, report = verify(work, mdir, "orthofinder_search.sh")
        assert ok is False
        assert "missing Blast1_1.txt.gz" in report
        assert "sbatch --array=1 orthofinder_search.sh" in report
        assert "did not run" in report


class TestMain:
    def test_complete_returns_zero(self, tmp_path):
        pairs = [(i, j) for i in range(2) for j in range(2)]
        work = _workdir(tmp_path, [0, 1], pairs)
        mdir = _manifests(tmp_path, {0: pairs})
        rc = main(["--work-dir", str(work), "--manifest-dir", str(mdir)])
        assert rc == 0

    def test_incomplete_returns_one(self, tmp_path):
        work = _workdir(tmp_path, [0, 1], [(0, 0), (0, 1), (1, 0)])
        mdir = _manifests(tmp_path, {0: [(0, 0), (0, 1)], 1: [(1, 0), (1, 1)]})
        rc = main(["--work-dir", str(work), "--manifest-dir", str(mdir)])
        assert rc == 1

    def test_missing_workdir_returns_one(self, tmp_path):
        rc = main(
            [
                "--work-dir", str(tmp_path / "nope"),
                "--manifest-dir", str(tmp_path / "m"),
            ]
        )
        assert rc == 1
