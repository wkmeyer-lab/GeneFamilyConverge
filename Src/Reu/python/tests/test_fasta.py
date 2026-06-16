"""
Tests for convgeno.io.fasta

Run with:  pytest tests/test_fasta.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.io.fasta import (
    FastaRecord,
    _parse_gene_id_ensembl,
    _parse_gene_id_ncbi,
    _strip_compression_suffix,
    detect_compression,
    detect_header_format,
    discover_fasta_files,
    filter_longest_isoforms,
    iter_fasta,
    open_fasta,
    parse_gene_id,
    process_directory,
)


# ===================================================================
#  1. detect_compression
# ===================================================================


class TestDetectCompression:
    def test_plain_text(self, ensembl_fasta_file: Path):
        assert detect_compression(ensembl_fasta_file) is None

    def test_gzip(self, gzip_fasta_file: Path):
        assert detect_compression(gzip_fasta_file) == "gzip"

    def test_bz2(self, bz2_fasta_file: Path):
        assert detect_compression(bz2_fasta_file) == "bz2"

    def test_gzip_without_extension(self, gzip_fasta_file: Path):
        """Magic bytes should detect gzip even without .gz suffix."""
        renamed = gzip_fasta_file.parent / "no_extension_file"
        gzip_fasta_file.rename(renamed)
        assert detect_compression(renamed) == "gzip"

    def test_nonexistent_file_with_gz_suffix(self, tmp_path: Path):
        """Suffix fallback when magic-byte read fails."""
        fake = tmp_path / "ghost.fa.gz"
        # Don't actually create the file — detect_compression catches OSError
        # and falls back to suffix.
        # However, the file doesn't exist so open() will raise.
        # detect_compression catches OSError, then checks suffix.
        # We need the file to exist for the open() to succeed.
        # Test suffix fallback with a non-gzip file that has .gz extension.
        fake.write_bytes(b"not actually gzip")
        assert detect_compression(fake) == "gzip"  # suffix fallback


# ===================================================================
#  2. open_fasta
# ===================================================================


class TestOpenFasta:
    def test_plain_text(self, ensembl_fasta_file: Path):
        with open_fasta(ensembl_fasta_file) as fh:
            first_line = next(iter(fh))
            assert first_line.startswith(">ENSP")

    def test_gzip(self, gzip_fasta_file: Path):
        with open_fasta(gzip_fasta_file) as fh:
            first_line = next(iter(fh))
            assert first_line.startswith(">ENSP")

    def test_bz2(self, bz2_fasta_file: Path):
        with open_fasta(bz2_fasta_file) as fh:
            first_line = next(iter(fh))
            assert first_line.startswith(">ENSP")

    def test_missing_file_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            open_fasta(tmp_path / "does_not_exist.fa")


# ===================================================================
#  3. iter_fasta (BioPython SeqIO wrapper)
# ===================================================================


class TestIterFasta:
    def test_ensembl_record_count(self, ensembl_fasta_file: Path):
        with open_fasta(ensembl_fasta_file) as fh:
            records = list(iter_fasta(fh))
        assert len(records) == 3

    def test_record_is_named_tuple(self, ensembl_fasta_file: Path):
        with open_fasta(ensembl_fasta_file) as fh:
            rec = next(iter_fasta(fh))
        assert isinstance(rec, FastaRecord)
        assert hasattr(rec, "header")
        assert hasattr(rec, "sequence")

    def test_header_has_no_angle_bracket(self, ensembl_fasta_file: Path):
        with open_fasta(ensembl_fasta_file) as fh:
            for rec in iter_fasta(fh):
                assert not rec.header.startswith(">")

    def test_sequence_has_no_whitespace(self, ensembl_fasta_file: Path):
        with open_fasta(ensembl_fasta_file) as fh:
            for rec in iter_fasta(fh):
                assert "\n" not in rec.sequence
                assert " " not in rec.sequence

    def test_gzip_streaming(self, gzip_fasta_file: Path):
        """Ensure SeqIO works through a gzip text handle."""
        with open_fasta(gzip_fasta_file) as fh:
            records = list(iter_fasta(fh))
        assert len(records) == 3

    def test_empty_file(self, tmp_path: Path):
        empty = tmp_path / "empty.fa"
        empty.write_text("", encoding="utf-8")
        with open_fasta(empty) as fh:
            records = list(iter_fasta(fh))
        assert records == []


# ===================================================================
#  4. Header format detection
# ===================================================================


class TestDetectHeaderFormat:
    def test_ensembl_gene_colon(self):
        h = "ENSP00000000001.1 pep gene:ENSG00000000001.1"
        assert detect_header_format(h) == "ensembl"

    def test_ensembl_gene_equals(self):
        h = "ENSP00000000001.1 pep gene=ENSG00000000001.1"
        assert detect_header_format(h) == "ensembl"

    def test_ncbi_xp(self):
        h = "XP_001000001.1 serine dehydratase [Mus musculus]"
        assert detect_header_format(h) == "ncbi"

    def test_ncbi_np(self):
        h = "NP_000001.1 alpha-1-B glycoprotein precursor [Homo sapiens]"
        assert detect_header_format(h) == "ncbi"

    def test_ncbi_without_brackets_falls_to_ensembl(self):
        # Missing [species] → not recognized as NCBI
        h = "XP_001000001.1 some protein"
        assert detect_header_format(h) == "ensembl"

    def test_ens_prefix_without_gene_tag(self):
        h = "ENSMUSP00000000001.1 pep transcript:ENSMUST00000000001.1"
        assert detect_header_format(h) == "ensembl"

    def test_ambiguous_defaults_to_ensembl(self):
        h = "some_random_accession a description"
        assert detect_header_format(h) == "ensembl"


# ===================================================================
#  5. Gene-ID parsing
# ===================================================================


class TestParseGeneIdEnsembl:
    def test_gene_colon(self):
        h = "ENSP00000000001.1 pep gene:ENSG00000000001.1 transcript:ENST01"
        assert _parse_gene_id_ensembl(h) == "ENSG00000000001.1"

    def test_gene_equals(self):
        h = "ENSP00000000001.1 pep gene=ENSG00000000001.1"
        assert _parse_gene_id_ensembl(h) == "ENSG00000000001.1"

    def test_no_gene_tag_returns_none(self):
        h = "ENSP00000000001.1 pep transcript:ENST00000000001.1"
        assert _parse_gene_id_ensembl(h) is None

    def test_locus_tag_not_matched(self):
        # gene: and gene= only — locus: should NOT match (different semantics)
        h = "ACC001 locus:LOC12345"
        assert _parse_gene_id_ensembl(h) is None


class TestParseGeneIdNcbi:
    def test_isoforms_collapse(self):
        h1 = "XP_001.1 serine dehydratase-like isoform X1 [Mus musculus]"
        h2 = "XP_002.1 serine dehydratase-like isoform X2 [Mus musculus]"
        assert _parse_gene_id_ncbi(h1) == _parse_gene_id_ncbi(h2)

    def test_paralogs_without_isoform_stay_distinct(self):
        h1 = "XP_003.1 40S ribosomal protein S12-like [Mus musculus]"
        h2 = "XP_004.1 40S ribosomal protein S12-like [Mus musculus]"
        assert _parse_gene_id_ncbi(h1) != _parse_gene_id_ncbi(h2)

    def test_isoform_with_letter_suffix(self):
        h1 = "XP_001.1 protein X isoform X1A [Species]"
        h2 = "XP_002.1 protein X isoform X2 [Species]"
        assert _parse_gene_id_ncbi(h1) == _parse_gene_id_ncbi(h2)


class TestParseGeneIdAuto:
    def test_auto_ensembl(self):
        h = "ENSP00000000001.1 pep gene:ENSG00000000001.1"
        assert parse_gene_id(h, "auto") == "ENSG00000000001.1"

    def test_auto_ncbi(self):
        h = "XP_001.1 serine dehydratase-like isoform X1 [Mus musculus]"
        key = parse_gene_id(h, "auto")
        assert key is not None
        assert "XP_001.1" not in key  # accession should be stripped

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError, match="Unknown header format"):
            parse_gene_id("anything", "bogus")  # type: ignore[arg-type]


# ===================================================================
#  6. filter_longest_isoforms
# ===================================================================


class TestFilterLongestIsoforms:
    def test_ensembl_basic(self, ensembl_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "out.fa"
        stats = filter_longest_isoforms(
            input_paths=[ensembl_fasta_file],
            output_path=out,
        )
        assert stats["total_records"] == 3
        assert stats["unique_genes"] == 2
        assert stats["unidentified"] == 0
        assert stats["output_records"] == 2
        assert stats["detected_format"] == "ensembl"
        assert out.exists()

        # Verify the longer isoform was kept.
        text = out.read_text()
        assert "MASEQKLISEEDLMASEQK" in text  # 19-char isoform
        assert text.count(">") == 2

    def test_ncbi_basic(self, ncbi_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "out.faa"
        stats = filter_longest_isoforms(
            input_paths=[ncbi_fasta_file],
            output_path=out,
        )
        # 2 isoforms collapse → 1 gene; 2 paralogs → 2 genes; total 3
        assert stats["unique_genes"] == 3
        assert stats["output_records"] == 3

    def test_gzip_input(self, gzip_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "out.fa"
        stats = filter_longest_isoforms(
            input_paths=[gzip_fasta_file], output_path=out
        )
        assert stats["unique_genes"] == 2

    def test_bz2_input(self, bz2_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "out.fa"
        stats = filter_longest_isoforms(
            input_paths=[bz2_fasta_file], output_path=out
        )
        assert stats["unique_genes"] == 2

    def test_multi_file_input(self, split_proteome_files: list[Path], tmp_path: Path):
        out = tmp_path / "merged.fa"
        stats = filter_longest_isoforms(
            input_paths=split_proteome_files, output_path=out
        )
        # chr1 has 1 record (gene1, len 13)
        # chr2 has 2 records (gene1 len 19, gene2 len 14)
        # After merging: gene1 → keep len 19, gene2 → keep len 14
        assert stats["total_records"] == 3
        assert stats["unique_genes"] == 2
        assert stats["output_records"] == 2

    def test_unidentified_records_preserved(
        self, ensembl_no_gene_file: Path, tmp_path: Path
    ):
        out = tmp_path / "out.fa"
        stats = filter_longest_isoforms(
            input_paths=[ensembl_no_gene_file], output_path=out
        )
        assert stats["unidentified"] == 1
        assert stats["output_records"] == 1
        # The record should still appear in output.
        assert out.read_text().count(">") == 1

    def test_original_headers_preserved(
        self, ensembl_fasta_file: Path, tmp_path: Path
    ):
        out = tmp_path / "out.fa"
        filter_longest_isoforms(
            input_paths=[ensembl_fasta_file], output_path=out
        )
        text = out.read_text()
        # The full original Ensembl header should be present.
        assert "gene:ENSG00000000001.1" in text
        assert "transcript:ENST00000000002.1" in text  # the longer isoform

    def test_output_sorted_by_gene_id(
        self, ensembl_fasta_file: Path, tmp_path: Path
    ):
        out = tmp_path / "out.fa"
        filter_longest_isoforms(
            input_paths=[ensembl_fasta_file], output_path=out
        )
        headers = [
            l for l in out.read_text().splitlines() if l.startswith(">")
        ]
        # ENSG00000000001 should come before ENSG00000000002
        assert "ENSG00000000001" in headers[0]
        assert "ENSG00000000002" in headers[1]

    def test_output_creates_parent_dirs(self, ensembl_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "deep" / "nested" / "out.fa"
        filter_longest_isoforms(
            input_paths=[ensembl_fasta_file], output_path=out
        )
        assert out.exists()


# ===================================================================
#  7. Duplicate accession handling
# ===================================================================


class TestDuplicateHandling:
    def test_error_mode_raises(self, duplicate_fasta_file: Path, tmp_path: Path):
        with pytest.raises(ValueError, match="Duplicate accession"):
            filter_longest_isoforms(
                input_paths=[duplicate_fasta_file],
                output_path=tmp_path / "out.fa",
                on_duplicate="error",
            )

    def test_warn_mode_continues(self, duplicate_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "out.fa"
        stats = filter_longest_isoforms(
            input_paths=[duplicate_fasta_file],
            output_path=out,
            on_duplicate="warn",
        )
        # Duplicate is dropped; first occurrence kept.
        assert stats["output_records"] == 1

    def test_skip_mode_continues(self, duplicate_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "out.fa"
        stats = filter_longest_isoforms(
            input_paths=[duplicate_fasta_file],
            output_path=out,
            on_duplicate="skip",
        )
        assert stats["output_records"] == 1


# ===================================================================
#  8. memory_mode='low' guard
# ===================================================================


class TestLowMemoryMode:
    """Verify that memory_mode='low' produces identical results to 'normal'."""

    def test_ensembl_matches_normal(self, ensembl_fasta_file: Path, tmp_path: Path):
        out_normal = tmp_path / "normal.fa"
        out_low = tmp_path / "low.fa"
        sn = filter_longest_isoforms(
            [ensembl_fasta_file], out_normal, memory_mode="normal",
        )
        sl = filter_longest_isoforms(
            [ensembl_fasta_file], out_low, memory_mode="low",
        )
        assert sn == sl
        assert out_normal.read_text() == out_low.read_text()

    def test_ncbi_matches_normal(self, ncbi_fasta_file: Path, tmp_path: Path):
        out_normal = tmp_path / "normal.faa"
        out_low = tmp_path / "low.faa"
        sn = filter_longest_isoforms(
            [ncbi_fasta_file], out_normal, memory_mode="normal",
        )
        sl = filter_longest_isoforms(
            [ncbi_fasta_file], out_low, memory_mode="low",
        )
        assert sn == sl
        assert out_normal.read_text() == out_low.read_text()

    def test_gzip_input(self, gzip_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "low_gz.fa"
        stats = filter_longest_isoforms(
            [gzip_fasta_file], out, memory_mode="low",
        )
        assert stats["unique_genes"] == 2
        assert out.exists()

    def test_bz2_input(self, bz2_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "low_bz2.fa"
        stats = filter_longest_isoforms(
            [bz2_fasta_file], out, memory_mode="low",
        )
        assert stats["unique_genes"] == 2

    def test_multi_file(self, split_proteome_files: list[Path], tmp_path: Path):
        out_normal = tmp_path / "normal.fa"
        out_low = tmp_path / "low.fa"
        sn = filter_longest_isoforms(
            split_proteome_files, out_normal, memory_mode="normal",
        )
        sl = filter_longest_isoforms(
            split_proteome_files, out_low, memory_mode="low",
        )
        assert sn == sl
        assert out_normal.read_text() == out_low.read_text()

    def test_duplicate_error(self, duplicate_fasta_file: Path, tmp_path: Path):
        with pytest.raises(ValueError, match="Duplicate accession"):
            filter_longest_isoforms(
                [duplicate_fasta_file],
                tmp_path / "x.fa",
                memory_mode="low",
                on_duplicate="error",
            )

    def test_duplicate_warn(self, duplicate_fasta_file: Path, tmp_path: Path):
        out = tmp_path / "low_warn.fa"
        stats = filter_longest_isoforms(
            [duplicate_fasta_file], out,
            memory_mode="low", on_duplicate="warn",
        )
        assert stats["output_records"] == 1

    def test_unidentified_preserved(
        self, ensembl_no_gene_file: Path, tmp_path: Path,
    ):
        out = tmp_path / "low_unid.fa"
        stats = filter_longest_isoforms(
            [ensembl_no_gene_file], out, memory_mode="low",
        )
        assert stats["unidentified"] == 1
        assert stats["output_records"] == 1


# ===================================================================
#  9. discover_fasta_files
# ===================================================================


class TestDiscoverFastaFiles:
    def test_finds_all_fasta(self, fasta_dir: Path):
        found = discover_fasta_files(fasta_dir)
        names = {f.name for f in found}
        assert "species_a.fa" in names
        assert "species_b.faa" in names
        assert "species_c.fasta.gz" in names

    def test_ignores_non_fasta(self, fasta_dir: Path):
        found = discover_fasta_files(fasta_dir)
        names = {f.name for f in found}
        assert "README.txt" not in names

    def test_returns_sorted(self, fasta_dir: Path):
        found = discover_fasta_files(fasta_dir)
        assert found == sorted(found)

    def test_not_a_directory_raises(self, ensembl_fasta_file: Path):
        with pytest.raises(NotADirectoryError):
            discover_fasta_files(ensembl_fasta_file)

    def test_empty_dir_returns_empty(self, empty_dir: Path):
        assert discover_fasta_files(empty_dir) == []


# ===================================================================
#  10. process_directory
# ===================================================================


class TestProcessDirectory:
    def test_batch_processing(self, fasta_dir: Path, tmp_path: Path):
        out_dir = tmp_path / "output"
        all_stats = process_directory(fasta_dir, out_dir)

        # Should have processed all three FASTA files.
        assert len(all_stats) == 3

        # Output files should exist with compression suffix stripped.
        assert (out_dir / "species_a.fa").exists()
        assert (out_dir / "species_b.faa").exists()
        assert (out_dir / "species_c.fasta").exists()

    def test_empty_dir_raises(self, empty_dir: Path, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="No FASTA files"):
            process_directory(empty_dir, tmp_path / "out")

    def test_creates_output_dir(self, fasta_dir: Path, tmp_path: Path):
        out_dir = tmp_path / "new" / "output"
        process_directory(fasta_dir, out_dir)
        assert out_dir.is_dir()


# ===================================================================
#  11. Helpers
# ===================================================================


class TestHelpers:
    @pytest.mark.parametrize(
        "name, expected",
        [
            ("species.fa.gz", "species.fa"),
            ("species.fasta.bz2", "species.fasta"),
            ("species.faa", "species.faa"),
            ("noext.gz", "noext"),
        ],
    )
    def test_strip_compression_suffix(self, name: str, expected: str):
        assert _strip_compression_suffix(name) == expected