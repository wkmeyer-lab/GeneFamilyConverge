"""
convgeno.io.paths — Cross-platform path resolution and directory helpers.

Functions to implement:
    ensure_dir(path)
        Create directory (and parents) if it doesn't exist.

    list_fastas(directory) -> list[Path]
        Return sorted list of .fasta / .fa files in a directory.

    species_name_from_filename(path) -> str
        Extract species name from a proteome FASTA filename.
        E.g., "Homo_sapiens.fasta" -> "Homo_sapiens"
"""

from pathlib import Path
