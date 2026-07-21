"""Tests for convgeno.io.tables."""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.io.tables import read_gene_counts, read_tsv, write_tsv


class TestReadWriteTsv:
    def test_round_trip(self, tmp_path: Path):
        rows = [
            {"a": "1", "b": "x"},
            {"a": "2", "b": "y"},
        ]
        p = tmp_path / "t.tsv"
        write_tsv(rows, p, ["a", "b"])
        assert read_tsv(p) == rows

    def test_uses_unix_line_endings(self, tmp_path: Path):
        p = tmp_path / "t.tsv"
        write_tsv([{"a": "1", "b": "2"}], p, ["a", "b"])
        assert "\r" not in p.read_text(encoding="utf-8")

    def test_extra_keys_ignored_and_missing_blank(self, tmp_path: Path):
        p = tmp_path / "t.tsv"
        write_tsv([{"a": "1", "z": "ignored"}], p, ["a", "b"])
        assert read_tsv(p) == [{"a": "1", "b": ""}]

    def test_read_missing_file_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            read_tsv(tmp_path / "nope.tsv")


class TestReadGeneCounts:
    def test_parses_species_and_counts(self, gene_count_file: Path):
        data = read_gene_counts(gene_count_file)
        assert data["species"] == ["human", "cat", "dog"]  # Total dropped
        assert data["families"] == ["OG0000000", "OG0000001", "OG0000002"]
        assert data["counts"]["OG0000000"] == {"human": 3, "cat": 2, "dog": 1}
        assert data["counts"]["OG0000002"]["human"] == 150

    def test_no_trailing_total_column(self, tmp_path: Path):
        p = tmp_path / "gc.tsv"
        p.write_text("Orthogroup\ta\tb\nOG0\t1\t2\n", encoding="utf-8")
        data = read_gene_counts(p)
        assert data["species"] == ["a", "b"]
        assert data["counts"]["OG0"] == {"a": 1, "b": 2}

    def test_non_integer_count_raises(self, tmp_path: Path):
        p = tmp_path / "bad.tsv"
        p.write_text("Orthogroup\ta\tb\tTotal\nOG0\t1\tNA\t1\n", encoding="utf-8")
        with pytest.raises(ValueError, match="Non-integer"):
            read_gene_counts(p)

    def test_missing_file_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            read_gene_counts(tmp_path / "nope.tsv")

    def test_empty_file_raises(self, tmp_path: Path):
        p = tmp_path / "empty.tsv"
        p.write_text("", encoding="utf-8")
        with pytest.raises(ValueError, match="empty"):
            read_gene_counts(p)
