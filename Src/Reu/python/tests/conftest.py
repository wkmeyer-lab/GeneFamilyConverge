"""
Shared pytest fixtures for convgeno tests.

All fixtures that create files use pytest's ``tmp_path`` so they are
automatically cleaned up after each test session.
"""

from __future__ import annotations

import bz2
import gzip
from pathlib import Path

import pytest

# ===================================================================
#  Raw FASTA content strings
# ===================================================================

ENSEMBL_FASTA = """\
>ENSP00000000001.1 pep chromosome:GRCh38:1:1000:2000:1 gene:ENSG00000000001.1 transcript:ENST00000000001.1
MASEQKLISEEDL
>ENSP00000000002.1 pep chromosome:GRCh38:1:1000:3000:1 gene:ENSG00000000001.1 transcript:ENST00000000002.1
MASEQKLISEEDLMASEQK
>ENSP00000000003.1 pep chromosome:GRCh38:2:5000:6000:1 gene:ENSG00000000002.1 transcript:ENST00000000003.1
MKTLLLTLVVVALA
"""

# Gene ENSG00000000001.1: two isoforms (lengths 13 and 19 → keep 19)
# Gene ENSG00000000002.1: one isoform  (length 14 → keep)
# Expected output: 2 genes, 3 total records, 2 written


NCBI_FASTA = """\
>XP_001000001.1 serine dehydratase-like isoform X1 [Mus musculus]
MASEQKLISEEDL
>XP_001000002.1 serine dehydratase-like isoform X2 [Mus musculus]
MASEQKLISEEDLMASEQK
>XP_001000003.1 40S ribosomal protein S12-like [Mus musculus]
MKTLLLTLVVVALA
>XP_001000004.1 40S ribosomal protein S12-like [Mus musculus]
MKTLLLTLVVVALARK
"""

# "serine dehydratase-like": two isoforms (isoform X1 / X2 → keep X2, len 19)
# "40S ribosomal protein S12-like": two paralogs (no isoform tag →
#     different accessions keep them distinct → both written)
# Expected output: 3 genes, 4 total records, 3 written


ENSEMBL_NO_GENE_TAG = """\
>SOMEACC001 pep chromosome:GRCh38:1:1000:2000:1 transcript:ENST00000000001.1
MASEQKLISEEDL
"""
# No gene: or gene= tag → unidentified


DUPLICATE_ACCESSION_FASTA = """\
>ENSP00000000001.1 pep gene:ENSG00000000001.1
MASEQKLISEEDL
>ENSP00000000001.1 pep gene:ENSG00000000001.1
MASEQKLISEEDLMORE
"""


# ===================================================================
#  Fixtures: plain-text FASTA files
# ===================================================================


@pytest.fixture()
def ensembl_fasta_file(tmp_path: Path) -> Path:
    """Write an Ensembl-format FASTA file and return its path."""
    p = tmp_path / "Homo_sapiens.pep.all.fa"
    p.write_text(ENSEMBL_FASTA, encoding="utf-8")
    return p


@pytest.fixture()
def ncbi_fasta_file(tmp_path: Path) -> Path:
    """Write an NCBI-format FASTA file and return its path."""
    p = tmp_path / "GCF_000001635.27_protein.faa"
    p.write_text(NCBI_FASTA, encoding="utf-8")
    return p


@pytest.fixture()
def ensembl_no_gene_file(tmp_path: Path) -> Path:
    p = tmp_path / "mystery.fasta"
    p.write_text(ENSEMBL_NO_GENE_TAG, encoding="utf-8")
    return p


@pytest.fixture()
def duplicate_fasta_file(tmp_path: Path) -> Path:
    p = tmp_path / "duplicated.fa"
    p.write_text(DUPLICATE_ACCESSION_FASTA, encoding="utf-8")
    return p


# ===================================================================
#  Fixtures: compressed FASTA files
# ===================================================================


@pytest.fixture()
def gzip_fasta_file(tmp_path: Path) -> Path:
    """Write a gzip-compressed Ensembl FASTA and return its path."""
    p = tmp_path / "Homo_sapiens.pep.all.fa.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(ENSEMBL_FASTA)
    return p


@pytest.fixture()
def bz2_fasta_file(tmp_path: Path) -> Path:
    """Write a bz2-compressed Ensembl FASTA and return its path."""
    p = tmp_path / "Homo_sapiens.pep.all.fa.bz2"
    with bz2.open(p, "wt", encoding="utf-8") as fh:
        fh.write(ENSEMBL_FASTA)
    return p


# ===================================================================
#  Fixtures: directories of FASTA files
# ===================================================================


@pytest.fixture()
def fasta_dir(tmp_path: Path) -> Path:
    """Create a directory with several FASTA files (mixed formats)."""
    d = tmp_path / "proteomes"
    d.mkdir()

    (d / "species_a.fa").write_text(ENSEMBL_FASTA, encoding="utf-8")
    (d / "species_b.faa").write_text(NCBI_FASTA, encoding="utf-8")

    with gzip.open(d / "species_c.fasta.gz", "wt", encoding="utf-8") as fh:
        fh.write(ENSEMBL_FASTA)

    # Non-FASTA file that should be ignored.
    (d / "README.txt").write_text("not a fasta file\n")

    return d


@pytest.fixture()
def empty_dir(tmp_path: Path) -> Path:
    d = tmp_path / "empty"
    d.mkdir()
    return d


# ===================================================================
#  Fixtures: split proteome (multi-file for one species)
# ===================================================================

ENSEMBL_SPLIT_CHR1 = """\
>ENSP00000000001.1 pep gene:ENSG00000000001.1 transcript:ENST00000000001.1
MASEQKLISEEDL
"""

ENSEMBL_SPLIT_CHR2 = """\
>ENSP00000000002.1 pep gene:ENSG00000000001.1 transcript:ENST00000000002.1
MASEQKLISEEDLMASEQK
>ENSP00000000003.1 pep gene:ENSG00000000002.1 transcript:ENST00000000003.1
MKTLLLTLVVVALA
"""


@pytest.fixture()
def split_proteome_files(tmp_path: Path) -> list[Path]:
    """Two FASTA files representing a split proteome for one species."""
    f1 = tmp_path / "species_chr1.fa"
    f2 = tmp_path / "species_chr2.fa"
    f1.write_text(ENSEMBL_SPLIT_CHR1, encoding="utf-8")
    f2.write_text(ENSEMBL_SPLIT_CHR2, encoding="utf-8")
    return [f1, f2]


# ===================================================================
#  Fixtures: species tree + alignment (for r8s / tree tests)
# ===================================================================

# A rooted species tree with FastTree-style internal support values
# (0.98, 1.0) that r8s cannot parse and must be stripped.
SPECIES_TREE_NEWICK = "((human:0.10,cat:0.12)0.98:0.05,dog:0.20)1.0:0.0;\n"

# Concatenated species-tree alignment: 3 records, each 12 columns wide
# (gaps included).  Its column count (12) is the value of ``nsites``.
SPECIES_TREE_ALIGNMENT = (
    ">human\nMKT-GAC--DEF\n>cat\nMKTAGAC--DEF\n>dog\nMKT-GACAADEF\n"
)


@pytest.fixture()
def species_tree_file(tmp_path: Path) -> Path:
    """Write a rooted species tree (with support values) and return its path."""
    p = tmp_path / "SpeciesTree_rooted.txt"
    p.write_text(SPECIES_TREE_NEWICK, encoding="utf-8")
    return p


@pytest.fixture()
def species_tree_alignment_file(tmp_path: Path) -> Path:
    """Write a 12-column concatenated alignment and return its path."""
    p = tmp_path / "SpeciesTreeAlignment.fa"
    p.write_text(SPECIES_TREE_ALIGNMENT, encoding="utf-8")
    return p


@pytest.fixture()
def orthofinder_output_dir(tmp_path: Path) -> Path:
    """Build a minimal OrthoFinder (MSA-mode) output tree and return its root.

    Layout::

        orthofinder_out/
          Results_test/
            Species_Tree/SpeciesTree_rooted.txt
            MultipleSequenceAlignments/SpeciesTreeAlignment.fa
    """
    root = tmp_path / "orthofinder_out"
    results = root / "Results_test"
    st_dir = results / "Species_Tree"
    st_dir.mkdir(parents=True)
    (st_dir / "SpeciesTree_rooted.txt").write_text(
        SPECIES_TREE_NEWICK, encoding="utf-8"
    )
    msa_dir = results / "MultipleSequenceAlignments"
    msa_dir.mkdir(parents=True)
    (msa_dir / "SpeciesTreeAlignment.fa").write_text(
        SPECIES_TREE_ALIGNMENT, encoding="utf-8"
    )
    og_dir = results / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups_GeneCount.tsv").write_text(
        GENE_COUNT_TSV, encoding="utf-8"
    )
    return root


# A dated/ultrametric tree: every root-to-tip path sums to 94.
ULTRAMETRIC_TREE_NEWICK = "((human:47,cat:47):47,dog:94);\n"

# OrthoFinder Orthogroups_GeneCount.tsv: species human/cat/dog + trailing Total.
# OG0000002 is a "large" family (>=100 in human) for size-filter tests.
GENE_COUNT_TSV = (
    "Orthogroup\thuman\tcat\tdog\tTotal\n"
    "OG0000000\t3\t2\t1\t6\n"
    "OG0000001\t0\t1\t2\t3\n"
    "OG0000002\t150\t1\t1\t152\n"
)


@pytest.fixture()
def ultrametric_tree_file(tmp_path: Path) -> Path:
    """Write a dated ultrametric tree and return its path."""
    p = tmp_path / "species_tree_ultrametric.nwk"
    p.write_text(ULTRAMETRIC_TREE_NEWICK, encoding="utf-8")
    return p


@pytest.fixture()
def gene_count_file(tmp_path: Path) -> Path:
    """Write an OrthoFinder-style gene count table and return its path."""
    p = tmp_path / "Orthogroups_GeneCount.tsv"
    p.write_text(GENE_COUNT_TSV, encoding="utf-8")
    return p
