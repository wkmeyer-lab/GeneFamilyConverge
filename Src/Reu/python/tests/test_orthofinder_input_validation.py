"""
Tests for convgeno.validation.orthofinder_inputs

Run with:  pytest tests/test_orthofinder_input_validation.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.validation.orthofinder_inputs import (
    ValidationResult,
    validate_orthofinder_inputs,
)


@pytest.fixture()
def make_fasta_dir(tmp_path: Path):
    """Factory fixture that creates a directory with fake FASTA files."""

    def _factory(
        species_names=None,
        extension=".fa",
        content=">gene1\nMKTLLIL\n",
        empty_files=None,
        bad_content_files=None,
    ):
        if species_names is None:
            species_names = [
                "Homo_sapiens", "Mus_musculus", "Bos_taurus",
                "Canis_familiaris", "Felis_catus",
            ]
        d = tmp_path / "proteomes"
        d.mkdir(exist_ok=True)
        empty_files = set(empty_files or [])
        bad_content_files = set(bad_content_files or [])
        for name in species_names:
            f = d / f"{name}{extension}"
            if name in empty_files:
                f.write_text("")
            elif name in bad_content_files:
                f.write_text("not a fasta file\njust plain text\n")
            else:
                f.write_text(content)
        return d

    return _factory


class TestValidDirectory:
    def test_passes(self, make_fasta_dir):
        d = make_fasta_dir()
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is True
        assert len(result.errors) == 0
        assert len(result.warnings) == 0
        assert len(result.fasta_files) == 5


class TestDirectoryProblems:
    def test_does_not_exist(self):
        result = validate_orthofinder_inputs("/nonexistent/path/proteomes")
        assert result.is_valid() is False
        assert any("does not exist" in e for e in result.errors)

    def test_is_a_file_not_directory(self, tmp_path: Path):
        f = tmp_path / "not_a_dir.fa"
        f.write_text(">seq\nACGT\n")
        result = validate_orthofinder_inputs(f)
        assert result.is_valid() is False
        assert any("not a directory" in e for e in result.errors)

    def test_no_fasta_files(self, tmp_path: Path):
        d = tmp_path / "empty_dir"
        d.mkdir()
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is False
        assert any("No FASTA files" in e for e in result.errors)


class TestSpeciesCounts:
    def test_too_few_species(self, make_fasta_dir):
        d = make_fasta_dir(species_names=["Sp1", "Sp2", "Sp3"])
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is False
        assert any("at least" in e for e in result.errors)

    def test_exactly_minimum_species(self, make_fasta_dir):
        d = make_fasta_dir(species_names=["Sp1", "Sp2", "Sp3", "Sp4"])
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is True


class TestFileContentChecks:
    def test_empty_fasta_file(self, make_fasta_dir):
        d = make_fasta_dir(
            species_names=["Sp1", "Sp2", "Sp3", "Sp4", "Empty_one"],
            empty_files=["Empty_one"],
        )
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is False
        assert any("empty" in e and "Empty_one" in e for e in result.errors)

    def test_file_not_starting_with_fasta_header(self, make_fasta_dir):
        d = make_fasta_dir(
            species_names=["Sp1", "Sp2", "Sp3", "Sp4", "Bad_one"],
            bad_content_files=["Bad_one"],
        )
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is False
        assert any("FASTA header" in e and "Bad_one" in e for e in result.errors)


class TestDuplicateSpecies:
    def test_duplicate_species_names(self, tmp_path: Path):
        d = tmp_path / "dup_test"
        d.mkdir()
        (d / "Homo_sapiens.fa").write_text(">g1\nMK\n")
        (d / "Homo_sapiens.fasta").write_text(">g2\nMK\n")
        (d / "Mus_musculus.fa").write_text(">g3\nMK\n")
        (d / "Bos_taurus.fa").write_text(">g4\nMK\n")
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is False
        assert any("Duplicate" in e and "Homo_sapiens" in e for e in result.errors)


class TestFilenameChecks:
    def test_unusual_characters_warning(self, make_fasta_dir):
        d = make_fasta_dir(
            species_names=["Sp1", "Sp2", "Sp3", "Sp 4 with spaces"]
        )
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is True
        assert len(result.warnings) > 0
        assert any("unusual characters" in w for w in result.warnings)


class TestExtensionHandling:
    def test_different_extensions_accepted(self, tmp_path: Path):
        d = tmp_path / "multi_ext"
        d.mkdir()
        (d / "Sp1.fa").write_text(">g1\nMK\n")
        (d / "Sp2.fasta").write_text(">g2\nMK\n")
        (d / "Sp3.faa").write_text(">g3\nMK\n")
        (d / "Sp4.fa").write_text(">g4\nMK\n")
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is True
        assert len(result.fasta_files) == 4

    def test_non_fasta_files_ignored(self, make_fasta_dir):
        d = make_fasta_dir()
        (d / "README.txt").write_text("notes")
        (d / "data.csv").write_text("a,b,c")
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is True
        assert len(result.fasta_files) == 5


class TestSummaryOutput:
    def test_summary_on_failure(self):
        result = validate_orthofinder_inputs("/nonexistent/path")
        s = result.summary()
        assert "VALIDATION FAILED" in s
        assert "ERROR" in s

    def test_summary_on_success(self, make_fasta_dir):
        d = make_fasta_dir()
        result = validate_orthofinder_inputs(d)
        s = result.summary()
        assert "Validation passed" in s
        assert "FASTA files found: 5" in s


class TestMultipleErrors:
    def test_multiple_errors_collected(self, tmp_path: Path):
        d = tmp_path / "multi_err"
        d.mkdir()
        (d / "Sp1.fa").write_text("")
        (d / "Sp2.fa").write_text("not fasta")
        result = validate_orthofinder_inputs(d)
        assert result.is_valid() is False
        assert len(result.errors) >= 3
