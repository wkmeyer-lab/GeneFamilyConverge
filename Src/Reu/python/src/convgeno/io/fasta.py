"""
convgeno.io.fasta
~~~~~~~~~~~~~~~~~

Reusable FASTA utilities for the convergent-evolution genomics pipeline.

Capabilities
------------
- Transparent decompression (gzip, bz2) via magic-byte detection
- Streaming FASTA record iteration backed by BioPython SeqIO
- Gene-ID extraction with auto-detection of Ensembl / NCBI header formats
- Single-pass longest-isoform filtering (v1; low-memory two-pass planned)
- Directory discovery of FASTA files by known bioinformatics extensions
- Batch processing of one-file-per-species directories

Dependencies
------------
- BioPython (``biopython >= 1.80``): used for robust FASTA parsing via
  ``Bio.SeqIO``.  The rest of the pipeline will also use NumPy (a
  BioPython transitive dependency), so this adds no extra weight.
"""

from __future__ import annotations

import bz2
import gzip
import logging
import re
import tempfile
from pathlib import Path
from typing import (
    Dict,
    Iterator,
    List,
    Literal,
    NamedTuple,
    Optional,
    TextIO,
)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

__all__ = [
    "open_fasta",
    "iter_fasta",
    "parse_gene_id",
    "detect_header_format",
    "filter_longest_isoforms",
    "discover_fasta_files",
    "process_directory",
    "FastaRecord",
    "HeaderFormat",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Base extensions recognised as FASTA.  ``.faa`` is the NCBI convention
#: for protein FASTA; ``.fa`` and ``.fasta`` are the generic conventions.
FASTA_EXTENSIONS: set[str] = {".fa", ".fasta", ".faa"}

#: Compression suffixes we handle.
COMPRESSED_SUFFIXES: set[str] = {".gz", ".bz2"}

#: FASTA sequence line width used when writing output (NCBI convention).
_LINE_WIDTH: int = 60

# Magic bytes ---------------------------------------------------------------
_GZIP_MAGIC = b"\x1f\x8b"
_BZ2_MAGIC = b"\x42\x5a\x68"  # "BZh"

# Header regexes ------------------------------------------------------------
#: Matches  gene:ENSG00000000001  or  gene=ENSG00000000001
_ENSEMBL_GENE_RE = re.compile(r"(?:gene[:=])(\S+)")

#: Matches  isoform X1 ,  isoform 2 ,  isoform X1A  etc.
_NCBI_ISOFORM_RE = re.compile(r"\s+isoform\s+X?\d+[A-Za-z]*")

# ---------------------------------------------------------------------------
# Public type aliases
# ---------------------------------------------------------------------------

HeaderFormat = Literal["ensembl", "ncbi", "auto"]


class FastaRecord(NamedTuple):
    """Lightweight record: original header (without ``>``) and sequence."""

    header: str
    sequence: str


# ===================================================================
#  1.  File opening with transparent decompression
# ===================================================================


def detect_compression(path: Path) -> Optional[str]:
    """Return ``'gzip'``, ``'bz2'``, or *None* for an uncompressed file.

    Detection strategy
    ------------------
    1. Read the first three bytes and compare against known magic bytes.
    2. If the magic bytes don't match, fall back to the file-name suffix
       (``.gz`` / ``.bz2``).

    This ordering means a file with no extension but valid gzip content
    is still detected correctly.
    """
    path = Path(path)
    try:
        with open(path, "rb") as fh:
            magic = fh.read(3)
    except OSError:
        magic = b""

    if len(magic) >= 2 and magic[:2] == _GZIP_MAGIC:
        return "gzip"
    if len(magic) >= 3 and magic[:3] == _BZ2_MAGIC:
        return "bz2"

    # Suffix fallback — handles the (rare) case where magic-byte read
    # failed but the name strongly suggests compression.
    if path.name.endswith(".gz"):
        return "gzip"
    if path.name.endswith(".bz2"):
        return "bz2"

    return None


def open_fasta(path: Path) -> TextIO:
    """Open a FASTA file and return a **text-mode** handle.

    Gzip and bz2 files are decompressed on-the-fly via the Python
    standard-library ``gzip`` / ``bz2`` modules, which stream
    line-by-line without loading the entire file into memory.

    Parameters
    ----------
    path : Path
        Path to a ``.fa``, ``.fasta``, or ``.faa`` file, optionally
        gzip- or bz2-compressed.

    Returns
    -------
    TextIO
        A file-like object yielding decoded text lines.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist on disk.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    compression = detect_compression(path)

    if compression == "gzip":
        return gzip.open(path, "rt", encoding="utf-8")
    if compression == "bz2":
        return bz2.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")  # noqa: SIM115


# ===================================================================
#  2.  Streaming FASTA parser (BioPython SeqIO)
# ===================================================================


def iter_fasta(handle: TextIO) -> Iterator[FastaRecord]:
    """Yield :class:`FastaRecord` tuples from an open text stream.

    Wraps :func:`Bio.SeqIO.parse` so that the rest of the codebase
    depends on a lightweight ``(header, sequence)`` tuple rather than
    the full ``SeqRecord`` object.  This keeps downstream code simple
    and makes it easy to swap the parser if ever needed.

    ``SeqIO.parse`` in FASTA mode is a generator—only one record is
    materialised at a time, so this is safe for arbitrarily large files
    when paired with a streaming handle from :func:`open_fasta`.
    """
    for record in SeqIO.parse(handle, "fasta"):
        # record.description is the full header line after '>',
        # including the accession and everything that follows.
        yield FastaRecord(
            header=record.description,
            sequence=str(record.seq),
        )


def iter_fasta_raw(handle: TextIO) -> Iterator[SeqRecord]:
    """Yield raw BioPython ``SeqRecord`` objects from a text stream.

    Use this when you need access to BioPython-specific features
    (e.g., ``record.seq`` as a ``Seq`` object for translation or
    alphabet checks).  For most pipeline uses, :func:`iter_fasta`
    returning lightweight :class:`FastaRecord` tuples is sufficient.
    """
    yield from SeqIO.parse(handle, "fasta")


# ===================================================================
#  3.  Gene-ID extraction (Ensembl / NCBI)
# ===================================================================


def detect_header_format(header: str) -> HeaderFormat:
    """Guess whether *header* follows Ensembl or NCBI conventions.

    Heuristics (evaluated in order)
    --------------------------------
    1. Contains ``gene:`` or ``gene=``  →  **ensembl**
    2. Accession starts with ``XP_``, ``NP_``, ``YP_``, or ``WP_``
       and the line ends with ``]`` (bracketed species name)  →  **ncbi**
    3. Accession starts with ``ENS``  →  **ensembl**
    4. Fallback  →  **ensembl** (the ``gene:`` tag is unambiguous, so
       trying Ensembl first is the safer default).
    """
    if "gene:" in header or "gene=" in header:
        return "ensembl"

    tokens = header.split()
    accession = tokens[0] if tokens else ""
    ncbi_prefixes = ("XP_", "NP_", "YP_", "WP_")
    if (
        any(accession.startswith(p) for p in ncbi_prefixes)
        and header.rstrip().endswith("]")
    ):
        return "ncbi"

    if accession.startswith("ENS"):
        return "ensembl"

    return "ensembl"


def _parse_gene_id_ensembl(header: str) -> Optional[str]:
    """Extract the gene identifier from an Ensembl-style header.

    Looks for a ``gene:ID`` or ``gene=ID`` token.  Returns *None* when
    no such token is present (the record will be treated as
    "unidentified" by the isoform filter).
    """
    m = _ENSEMBL_GENE_RE.search(header)
    return m.group(1) if m else None


def _parse_gene_id_ncbi(header: str) -> str:
    """Derive a gene-identity key from an NCBI-style header.

    Rules
    -----
    - Records **with** an ``isoform …`` tag are considered isoforms of
      the same gene.  The accession is stripped so that different
      accessions sharing the same description collapse to one key.
    - Records **without** an ``isoform`` tag retain their accession,
      so that paralogs with identical functional descriptions (e.g.,
      two distinct "40S ribosomal protein S12-like" genes) stay
      separate.

    This matches the logic used by OrthoFinder's
    ``primary_transcript.py``.
    """
    has_isoform = _NCBI_ISOFORM_RE.search(header) is not None
    cleaned = _NCBI_ISOFORM_RE.sub("", header)

    if has_isoform:
        # Drop the accession so all isoforms collapse to one key.
        parts = cleaned.split(None, 1)
        return parts[1] if len(parts) > 1 else cleaned
    # Keep accession → distinct paralogs stay separate.
    return cleaned


def parse_gene_id(header: str, fmt: HeaderFormat = "auto") -> Optional[str]:
    """Return the gene identifier for a FASTA header.

    Parameters
    ----------
    header : str
        Header text **without** the leading ``>``.
    fmt : HeaderFormat
        ``'ensembl'``, ``'ncbi'``, or ``'auto'`` (detect from the
        header itself via :func:`detect_header_format`).

    Returns
    -------
    str or None
        The gene identity key, or *None* if the gene could not be
        identified (Ensembl header with no ``gene:`` tag).
    """
    if fmt == "auto":
        fmt = detect_header_format(header)

    if fmt == "ensembl":
        return _parse_gene_id_ensembl(header)
    if fmt == "ncbi":
        return _parse_gene_id_ncbi(header)

    raise ValueError(f"Unknown header format: {fmt!r}")


# ===================================================================
#  4.  Longest-isoform filter
# ===================================================================


def _write_fasta_record(fh: TextIO, header: str, sequence: str) -> None:
    """Write one FASTA record with sequence wrapped at ``_LINE_WIDTH``."""
    fh.write(f">{header}\n")
    for i in range(0, len(sequence), _LINE_WIDTH):
        fh.write(sequence[i : i + _LINE_WIDTH])
        fh.write("\n")


def _read_record_at(fh: TextIO, offset: int) -> tuple[str, str]:
    """Seek to *offset* in a text-mode handle and read one FASTA record.

    Returns ``(header, sequence)`` where *header* has no leading ``>``.
    """
    fh.seek(offset)
    header_line = fh.readline()
    if not header_line or not header_line.startswith(">"):
        raise ValueError(
            f"Expected '>' at offset {offset}, got: {header_line!r}"
        )
    header = header_line.rstrip("\n\r")[1:]
    seq_parts: list[str] = []
    while True:
        peek = fh.tell()
        line = fh.readline()
        if not line or line.startswith(">"):
            fh.seek(peek)
            break
        seq_parts.append(line.rstrip("\n\r"))
    return header, "".join(seq_parts)


def _filter_low_memory(
    input_paths: List[Path],
    output_path: Path,
    header_format: HeaderFormat,
    on_duplicate: Literal["error", "warn", "skip"],
) -> dict:
    """Two-pass low-memory longest-isoform filter.

    Compressed inputs are decompressed to a temporary directory
    (gzip / bz2 streams are not seekable).  Pass 1 records only
    ``gene_id → (length, file_index, byte_offset)``; no sequences
    are held in memory.  Pass 2 seeks to winning offsets and writes.
    """
    with tempfile.TemporaryDirectory(prefix="convgeno_") as tmp_str:
        tmp_dir = Path(tmp_str)

        # --- Prepare seekable plain-text copies ---
        seekable: list[Path] = []
        for i, p in enumerate(input_paths):
            p = Path(p)
            if detect_compression(p) is not None:
                decomp = tmp_dir / f"decomp_{i}.fa"
                logger.info("Decompressing %s → %s", p.name, decomp.name)
                with open_fasta(p) as src, \
                     open(decomp, "w", encoding="utf-8") as dst:
                    for line in src:
                        dst.write(line)
                seekable.append(decomp)
            else:
                seekable.append(p)

        # --- PASS 1: scan gene IDs, lengths, file offsets ---
        best: Dict[str, tuple[int, int, int]] = {}  # gene → (len, fidx, off)
        seen_accessions: set[str] = set()
        unidentified_offsets: list[tuple[int, int]] = []

        total_records = 0
        unidentified_count = 0
        resolved_fmt: Optional[HeaderFormat] = (
            None if header_format == "auto" else header_format
        )

        for file_idx, fpath in enumerate(seekable):
            logger.info("Pass 1 — scanning %s", fpath.name)
            with open(fpath, "r", encoding="utf-8") as fh:
                while True:
                    record_offset = fh.tell()
                    line = fh.readline()
                    if not line:
                        break
                    if not line.startswith(">"):
                        continue

                    total_records += 1
                    header = line.rstrip("\n\r")[1:]

                    # --- duplicate check ---
                    accession = header.split()[0]
                    if accession in seen_accessions:
                        msg = (
                            f"Duplicate accession '{accession}' "
                            f"in {fpath.name}"
                        )
                        if on_duplicate == "error":
                            raise ValueError(msg)
                        if on_duplicate == "warn":
                            logger.warning(msg)
                        # Advance past sequence lines.
                        while True:
                            peek = fh.tell()
                            sl = fh.readline()
                            if not sl or sl.startswith(">"):
                                fh.seek(peek)
                                break
                        continue
                    seen_accessions.add(accession)

                    # --- auto-detect format from first record ---
                    if resolved_fmt is None:
                        resolved_fmt = detect_header_format(header)
                        logger.info(
                            "Auto-detected header format: %s", resolved_fmt
                        )

                    gene_id = parse_gene_id(header, resolved_fmt)

                    # Count sequence length WITHOUT storing sequence.
                    seq_len = 0
                    while True:
                        peek = fh.tell()
                        sl = fh.readline()
                        if not sl or sl.startswith(">"):
                            fh.seek(peek)
                            break
                        seq_len += len(sl.rstrip("\n\r"))

                    if gene_id is None:
                        unidentified_count += 1
                        unidentified_offsets.append(
                            (file_idx, record_offset)
                        )
                    elif gene_id not in best or seq_len > best[gene_id][0]:
                        best[gene_id] = (
                            seq_len, file_idx, record_offset,
                        )

        # --- PASS 2: seek to winners, read, write ---
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        handles: list[TextIO] = []
        try:
            handles = [
                open(p, "r", encoding="utf-8") for p in seekable
            ]
            output_records = 0
            with open(output_path, "w", encoding="utf-8") as out:
                for gene_id in sorted(best):
                    _, fidx, off = best[gene_id]
                    hdr, seq = _read_record_at(handles[fidx], off)
                    _write_fasta_record(out, hdr, seq)
                    output_records += 1

                for fidx, off in unidentified_offsets:
                    hdr, seq = _read_record_at(handles[fidx], off)
                    _write_fasta_record(out, hdr, seq)
                    output_records += 1
        finally:
            for h in handles:
                h.close()

    logger.info(
        "Low-memory: %d records -> %d genes + %d unidentified = %d written",
        total_records, len(best), unidentified_count, output_records,
    )

    return {
        "total_records": total_records,
        "unique_genes": len(best),
        "unidentified": unidentified_count,
        "output_records": output_records,
        "detected_format": resolved_fmt,
    }


def filter_longest_isoforms(
    input_paths: List[Path],
    output_path: Path,
    header_format: HeaderFormat = "auto",
    on_duplicate: Literal["error", "warn", "skip"] = "error",
    memory_mode: Literal["normal", "low"] = "normal",
) -> dict:
    """Keep only the longest isoform per gene across one or more FASTA files.

    Algorithm (``memory_mode='normal'``, the v1 default)
    ----------------------------------------------------
    Single-pass, dictionary-based.  For every gene encountered, the
    header + sequence of the longest isoform seen so far is held in a
    dict.  Peak memory ≈ total size of the filtered output (typically
    10–15 MB for a mammalian proteome).  The file is read exactly once.

    ``memory_mode='low'`` is for a two-pass strategy
    that records only file offsets on the first pass.  It will require
    decompressing twice for gzip/bz2 inputs. 

    Parameters
    ----------
    input_paths : list[Path]
        One or more FASTA files for a **single species**.  Multiple
        files are processed as a concatenated stream (useful when a
        proteome is split across per-chromosome downloads).
    output_path : Path
        Where to write the filtered FASTA.
    header_format : HeaderFormat
        ``'ensembl'``, ``'ncbi'``, or ``'auto'``.  When ``'auto'``,
        the format is detected once from the first record's header and
        reused for all subsequent records.
    on_duplicate : ``'error'`` | ``'warn'`` | ``'skip'``
        Policy for duplicate accessions (first whitespace-delimited
        token).  ``'error'`` (the default) raises ``ValueError``.
    memory_mode : ``'normal'`` | ``'low'``
        ``'normal'``: single-pass, dict in memory.
        ``'low'``: two-pass, offset-based.  Compressed inputs are
        decompressed to a temp directory first (gzip/bz2 streams
        are not seekable).  Sequences are never held in memory.

    Returns
    -------
    dict
        ``{ "total_records", "unique_genes", "unidentified",
           "output_records", "detected_format" }``
    """
    if memory_mode == "low":
        return _filter_low_memory(
            input_paths, output_path, header_format, on_duplicate,
        )

    # gene_id → (seq_length, original_header, sequence)
    best: Dict[str, tuple[int, str, str]] = {}
    seen_accessions: set[str] = set()

    # Unidentified records (no gene ID extracted) are preserved as-is.
    unidentified_records: list[tuple[str, str]] = []

    total_records = 0
    unidentified_count = 0

    # Resolve format once from the first header, not per-record.
    resolved_fmt: Optional[HeaderFormat] = (
        None if header_format == "auto" else header_format
    )

    for fasta_path in input_paths:
        fasta_path = Path(fasta_path)
        logger.info("Reading %s", fasta_path)

        with open_fasta(fasta_path) as handle:
            for record in iter_fasta(handle):
                total_records += 1

                # --- duplicate-accession check ---
                accession = record.header.split()[0]
                if accession in seen_accessions:
                    msg = (
                        f"Duplicate accession '{accession}' "
                        f"in {fasta_path.name}"
                    )
                    if on_duplicate == "error":
                        raise ValueError(msg)
                    if on_duplicate == "warn":
                        logger.warning(msg)
                    # Both 'warn' and 'skip' drop the duplicate.
                    continue
                seen_accessions.add(accession)

                # --- auto-detect format from first record ---
                if resolved_fmt is None:
                    resolved_fmt = detect_header_format(record.header)
                    logger.info(
                        "Auto-detected header format: %s", resolved_fmt
                    )

                # --- extract gene ID ---
                gene_id = parse_gene_id(record.header, resolved_fmt)

                if gene_id is None:
                    unidentified_count += 1
                    unidentified_records.append(
                        (record.header, record.sequence)
                    )
                    continue

                # --- keep longest ---
                seq_len = len(record.sequence)
                if gene_id not in best or seq_len > best[gene_id][0]:
                    best[gene_id] = (seq_len, record.header, record.sequence)

    # ---- write output ----
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    output_records = 0
    with open(output_path, "w", encoding="utf-8") as out:
        # Identified genes — sorted by gene ID for reproducible output.
        for gene_id in sorted(best):
            _, header, sequence = best[gene_id]
            _write_fasta_record(out, header, sequence)
            output_records += 1

        # Unidentified records — preserved with original headers.
        for header, sequence in unidentified_records:
            _write_fasta_record(out, header, sequence)
            output_records += 1

    logger.info(
        "Processed %d records -> %d genes + %d unidentified = %d written",
        total_records,
        len(best),
        unidentified_count,
        output_records,
    )

    return {
        "total_records": total_records,
        "unique_genes": len(best),
        "unidentified": unidentified_count,
        "output_records": output_records,
        "detected_format": resolved_fmt,
    }


# ===================================================================
#  5.  Directory discovery & batch processing
# ===================================================================


def discover_fasta_files(directory: Path) -> list[Path]:
    """Return a sorted list of FASTA files in *directory* (non-recursive).

    Recognised extensions: ``.fa``, ``.fasta``, ``.faa``, each
    optionally followed by ``.gz`` or ``.bz2``.

    Files that do not match these extensions are silently ignored.
    Rename ``.txt`` or other non-standard FASTA files before running.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise NotADirectoryError(f"Not a directory: {directory}")

    found: list[Path] = []
    for child in directory.iterdir():
        if child.is_dir():
            continue
        # Strip a trailing compression suffix to inspect the base ext.
        stem = child.name
        for csuf in COMPRESSED_SUFFIXES:
            if stem.endswith(csuf):
                stem = stem[: -len(csuf)]
                break
        if Path(stem).suffix.lower() in FASTA_EXTENSIONS:
            found.append(child)

    return sorted(found)


def _strip_compression_suffix(name: str) -> str:
    """Remove ``.gz`` or ``.bz2`` from the end of a filename."""
    for csuf in COMPRESSED_SUFFIXES:
        if name.endswith(csuf):
            return name[: -len(csuf)]
    return name


def process_directory(
    input_dir: Path,
    output_dir: Path,
    header_format: HeaderFormat = "auto",
    on_duplicate: Literal["error", "warn", "skip"] = "error",
) -> dict[str, dict]:
    """Filter longest isoforms for every FASTA file in *input_dir*.

    Each file is assumed to contain the proteome of **one species**.
    Output files keep the original filename (with any compression
    suffix removed) and are written into *output_dir*.

    This function processes files **sequentially**.  Cross-species
    parallelism is the workflow manager's job (Snakemake rule or SLURM
    array).

    Returns
    -------
    dict
        ``{ input_filename: stats_dict }`` for every processed file.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files found in {input_dir}. "
            f"Expected extensions: {', '.join(sorted(FASTA_EXTENSIONS))} "
            f"(optionally compressed with .gz or .bz2)."
        )

    logger.info("Found %d FASTA files in %s", len(fasta_files), input_dir)

    all_stats: dict[str, dict] = {}
    for fasta_path in fasta_files:
        out_name = _strip_compression_suffix(fasta_path.name)
        stats = filter_longest_isoforms(
            input_paths=[fasta_path],
            output_path=output_dir / out_name,
            header_format=header_format,
            on_duplicate=on_duplicate,
        )
        all_stats[fasta_path.name] = stats

    return all_stats