"""File I/O utilities for FASTA, tables, and path management."""

from convgeno.io.fasta import (
    discover_fasta_files,
    filter_longest_isoforms,
    iter_fasta,
    open_fasta,
    parse_gene_id,
    process_directory,
)

__all__ = [
    "open_fasta",
    "iter_fasta",
    "parse_gene_id",
    "filter_longest_isoforms",
    "discover_fasta_files",
    "process_directory",
]