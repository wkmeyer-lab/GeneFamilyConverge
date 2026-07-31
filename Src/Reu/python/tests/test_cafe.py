"""Tests for convgeno.external.cafe."""

from __future__ import annotations

from pathlib import Path

from convgeno.external.cafe import (
    build_command,
    format_gene_counts_for_cafe,
    validate_input,
    validate_lambda_tree,
    validate_output,
)
from convgeno.io.tables import read_tsv


class TestFormatGeneCountsForCafe:
    def test_reformats_header_and_counts(self, gene_count_file: Path, tmp_path: Path):
        out = tmp_path / "cafe_input.tsv"
        stats = format_gene_counts_for_cafe(gene_count_file, out)

        header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert header == ["Desc", "Family ID", "human", "cat", "dog"]  # Total dropped

        rows = read_tsv(out)
        assert stats["n_families"] == 3
        assert stats["n_written"] == 3
        assert stats["species"] == ["human", "cat", "dog"]
        assert rows[0]["Desc"] == "(null)"
        assert rows[0]["Family ID"] == "OG0000000"
        assert rows[0]["human"] == "3"

    def test_custom_description(self, gene_count_file: Path, tmp_path: Path):
        out = tmp_path / "c.tsv"
        format_gene_counts_for_cafe(gene_count_file, out, description="fam")
        assert read_tsv(out)[0]["Desc"] == "fam"

    def test_size_filter_excludes_large_families(
        self, gene_count_file: Path, tmp_path: Path
    ):
        out = tmp_path / "cafe_input.tsv"
        stats = format_gene_counts_for_cafe(gene_count_file, out, max_family_size=100)

        assert stats["n_written"] == 2  # OG0000002 (150) excluded
        assert stats["n_large_excluded"] == 1
        large = Path(stats["large_family_file"])
        assert large.is_file()
        assert read_tsv(large)[0]["Family ID"] == "OG0000002"


class TestBuildCommand:
    def test_minimal(self):
        assert build_command("counts.tsv", "tree.nwk") == [
            "cafe5",
            "-i",
            "counts.tsv",
            "-t",
            "tree.nwk",
        ]

    def test_with_gamma_and_extra_args(self):
        cmd = build_command(
            "c.tsv",
            "t.nwk",
            n_gamma_cats=5,
            extra_args="-p -o out",
            tool_path="/x/cafe5",
        )
        assert cmd[0] == "/x/cafe5"
        assert "-k" in cmd and "5" in cmd
        assert cmd[-3:] == ["-p", "-o", "out"]

    def test_with_lambda_tree(self):
        cmd = build_command("c.tsv", "t.nwk", lambda_tree="lambda.nwk")
        assert cmd == ["cafe5", "-i", "c.tsv", "-t", "t.nwk", "-y", "lambda.nwk"]

    def test_lambda_tree_precedes_gamma(self):
        cmd = build_command("c.tsv", "t.nwk", n_gamma_cats=3, lambda_tree="l.nwk")
        assert cmd == ["cafe5", "-i", "c.tsv", "-t", "t.nwk", "-y", "l.nwk", "-k", "3"]


class TestValidateInput:
    def _cafe_counts(self, gene_count_file: Path, tmp_path: Path) -> Path:
        out = tmp_path / "cafe_input.tsv"
        format_gene_counts_for_cafe(gene_count_file, out)
        return out

    def test_valid_pair(self, gene_count_file, ultrametric_tree_file, tmp_path):
        counts = self._cafe_counts(gene_count_file, tmp_path)
        assert validate_input(counts, ultrametric_tree_file) == []

    def test_non_ultrametric_tree_reported(
        self, gene_count_file, species_tree_file, tmp_path
    ):
        counts = self._cafe_counts(gene_count_file, tmp_path)
        errors = validate_input(counts, species_tree_file)
        assert any("ultrametric" in e.lower() for e in errors)

    def test_species_mismatch_reported(self, gene_count_file, tmp_path):
        counts = self._cafe_counts(gene_count_file, tmp_path)
        # tree tips: human, cat, frog -> dog missing, frog extra
        tree = tmp_path / "mismatch.nwk"
        tree.write_text("((human:47,cat:47):47,frog:94);\n", encoding="utf-8")
        errors = validate_input(counts, tree)
        joined = " ".join(errors)
        assert "frog" in joined and "dog" in joined

    def test_missing_count_file_reported(self, ultrametric_tree_file, tmp_path):
        errors = validate_input(tmp_path / "nope.tsv", ultrametric_tree_file)
        assert any("not found" in e.lower() for e in errors)


class TestValidateOutput:
    def test_missing_dir(self, tmp_path: Path):
        errors = validate_output(tmp_path / "nope")
        assert errors and "not found" in errors[0].lower()

    def test_no_results_file(self, tmp_path: Path):
        d = tmp_path / "out"
        d.mkdir()
        errors = validate_output(d)
        assert any("results.txt" in e for e in errors)

    def test_results_present(self, tmp_path: Path):
        d = tmp_path / "out"
        d.mkdir()
        (d / "Base_results.txt").write_text("ok", encoding="utf-8")
        assert validate_output(d) == []


class TestValidateLambdaTree:
    def _w(self, tmp_path: Path, name: str, text: str) -> Path:
        p = tmp_path / name
        p.write_text(text, encoding="utf-8")
        return p

    def test_valid(self, tmp_path: Path):
        t = self._w(tmp_path, "t.nwk", "((human:47,cat:47):47,dog:94);\n")
        y = self._w(tmp_path, "y.nwk", "((human:1,cat:1):1,dog:2);\n")
        assert validate_lambda_tree(y, t) == []

    def test_tip_mismatch_reported(self, tmp_path: Path):
        t = self._w(tmp_path, "t.nwk", "((human:47,cat:47):47,dog:94);\n")
        y = self._w(tmp_path, "y.nwk", "((human:1,cat:1):1,frog:2);\n")
        joined = " ".join(validate_lambda_tree(y, t))
        assert "frog" in joined and "dog" in joined

    def test_non_integer_labels_reported(self, tmp_path: Path):
        t = self._w(tmp_path, "t.nwk", "((human:47,cat:47):47,dog:94);\n")
        y = self._w(tmp_path, "y.nwk", "((human:1.5,cat:1):1,dog:2);\n")
        assert any("integer" in e.lower() for e in validate_lambda_tree(y, t))
