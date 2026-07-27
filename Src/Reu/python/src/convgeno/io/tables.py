"""
convgeno.io.tables
~~~~~~~~~~~~~~~~~~

Read and write TSV tables consistently, and parse OrthoFinder's gene-count
table.  Thin wrappers around the standard-library :mod:`csv` module — pandas is
deliberately kept out of the core library so ``convgeno`` stays lightweight
(Loc scripts may use pandas if they wish).

All output uses ``\\n`` line terminators (not ``\\r\\n``) so downstream tools
(CAFE-5, R) read the files cleanly on every platform.
"""

from __future__ import annotations

import csv
from pathlib import Path

__all__ = ["read_tsv", "write_tsv", "read_gene_counts"]


def read_tsv(path: Path | str) -> list[dict]:
    """Read a tab-delimited file into a list of row dicts (header-keyed)."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"TSV file not found: {path}")
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [dict(row) for row in reader]


def write_tsv(rows: list[dict], path: Path | str, columns: list[str]) -> None:
    """Write *rows* to a tab-delimited file with the given *columns* order.

    A header line is written first.  Keys in a row that are not in *columns*
    are ignored; missing keys are written as empty strings.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_gene_counts(path: Path | str) -> dict:
    """Parse OrthoFinder's ``Orthogroups_GeneCount.tsv`` into a structured dict.

    The file's header is ``Orthogroup<TAB>species1<TAB>...<TAB>Total``; the
    trailing ``Total`` column (if present) is dropped.

    Returns
    -------
    dict
        ``{"family_column": str, "species": list[str], "families": list[str],
        "counts": {family_id: {species: int}}}``.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file is empty, malformed, or has a non-integer count.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Gene count file not found: {path}")

    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Gene count file is empty: {path}") from exc

        if len(header) < 2:
            raise ValueError(
                f"Gene count file needs at least an ID column and one species "
                f"column: {path}"
            )

        family_column = header[0]
        species = header[1:]
        if species and species[-1].strip().lower() == "total":
            species = species[:-1]
        if not species:
            raise ValueError(f"No species columns found in {path}")

        families: list[str] = []
        counts: dict[str, dict[str, int]] = {}
        for lineno, row in enumerate(reader, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            family = row[0]
            fam_counts: dict[str, int] = {}
            for i, sp in enumerate(species):
                try:
                    fam_counts[sp] = int(row[i + 1])
                except (IndexError, ValueError) as exc:
                    raise ValueError(
                        f"Non-integer or missing count at line {lineno}, "
                        f"species '{sp}' in {path}"
                    ) from exc
            families.append(family)
            counts[family] = fam_counts

    return {
        "family_column": family_column,
        "species": species,
        "families": families,
        "counts": counts,
    }
