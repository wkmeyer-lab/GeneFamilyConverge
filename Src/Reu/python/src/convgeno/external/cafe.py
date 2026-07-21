"""
convgeno.external.cafe
~~~~~~~~~~~~~~~~~~~~~~~

Prepare inputs for and wrap the CAFE-5 CLI.

CAFE-5 needs two inputs:

1. A **tab-delimited count file** whose header is
   ``Desc<TAB>Family ID<TAB><species...>`` and whose rows hold a description,
   a family ID, and non-negative integer gene counts per species.
2. A **rooted, binary, ultrametric** Newick species tree whose tips match the
   count file's species columns.

This module reformats OrthoFinder's ``Orthogroups_GeneCount.tsv`` into that
layout, validates the pair of inputs, and builds the CAFE-5 command.  Actual
execution goes through :func:`convgeno.utils.command_runner.run`.
"""

from __future__ import annotations

import logging
import shlex
from pathlib import Path

from convgeno.io.tables import read_gene_counts, read_tsv, write_tsv
from convgeno.validation.trees import (
    check_tips_match_species,
    is_ultrametric,
    validate_tree,
)

logger = logging.getLogger(__name__)

__all__ = [
    "format_gene_counts_for_cafe",
    "validate_input",
    "validate_output",
    "build_command",
]

#: The two leading columns CAFE-5 expects before the per-species counts.
CAFE_LEADING_COLUMNS = ["Desc", "Family ID"]


def format_gene_counts_for_cafe(
    orthofinder_counts_path: Path | str,
    output_path: Path | str,
    *,
    description: str = "(null)",
    max_family_size: int | None = None,
) -> dict:
    """Reformat OrthoFinder gene counts into a CAFE-5 count file.

    Drops the trailing ``Total`` column, prepends the ``Desc`` and ``Family ID``
    columns, and writes a tab-delimited table.

    Parameters
    ----------
    orthofinder_counts_path : Path or str
        OrthoFinder ``Orthogroups_GeneCount.tsv``.
    output_path : Path or str
        Where to write the CAFE-5 count file.
    description : str
        Value for the ``Desc`` column (CAFE-5 uses ``(null)`` when unknown).
    max_family_size : int, optional
        If set, families in which any species has a count ``>=`` this value are
        excluded (mirrors the CAFE-5 tutorial's size filter) and written to a
        sibling ``*.large.tsv`` file.

    Returns
    -------
    dict
        Summary stats: ``n_families``, ``n_written``, ``n_species``,
        ``species`` (and ``n_large_excluded`` / ``large_family_file`` when
        filtering).
    """
    data = read_gene_counts(orthofinder_counts_path)
    species = data["species"]
    columns = [*CAFE_LEADING_COLUMNS, *species]

    kept: list[dict] = []
    large: list[dict] = []
    for family in data["families"]:
        fam_counts = data["counts"][family]
        row = {"Desc": description, "Family ID": family}
        row.update({sp: fam_counts[sp] for sp in species})
        if max_family_size is not None and any(
            fam_counts[sp] >= max_family_size for sp in species
        ):
            large.append(row)
        else:
            kept.append(row)

    output_path = Path(output_path)
    write_tsv(kept, output_path, columns)
    logger.info("Wrote %d gene families to %s", len(kept), output_path)

    stats: dict = {
        "n_families": len(data["families"]),
        "n_written": len(kept),
        "n_species": len(species),
        "species": list(species),
    }
    if max_family_size is not None:
        large_path = output_path.with_name(
            f"{output_path.stem}.large{output_path.suffix}"
        )
        write_tsv(large, large_path, columns)
        stats["n_large_excluded"] = len(large)
        stats["large_family_file"] = str(large_path)
        logger.info(
            "Excluded %d large families (>=%d copies) to %s",
            len(large),
            max_family_size,
            large_path,
        )
    return stats


def validate_input(count_file: Path | str, tree_file: Path | str) -> list[str]:
    """Validate a CAFE-5 count file + species tree pair.

    Checks the tree is parseable/rooted/binary and ultrametric, the count file
    has the expected header and non-negative integer counts, and the tree tips
    exactly match the count file's species columns.

    Returns
    -------
    list[str]
        Error messages; empty means the inputs are valid.
    """
    errors: list[str] = []

    errors.extend(validate_tree(tree_file))
    if not is_ultrametric(tree_file):
        errors.append(
            "Species tree is not ultrametric; CAFE-5 requires an ultrametric tree."
        )

    count_file = Path(count_file)
    if not count_file.is_file():
        errors.append(f"Count file not found: {count_file}")
        return errors

    with open(count_file, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if len(header) < 3 or header[0] != "Desc" or header[1] != "Family ID":
        errors.append(
            "Count file header must be 'Desc<TAB>Family ID<TAB>' followed by "
            "one column per species."
        )
        return errors

    species = header[2:]
    for i, row in enumerate(read_tsv(count_file), start=2):
        for sp in species:
            value = row.get(sp, "")
            try:
                count = int(value)
            except ValueError:
                errors.append(
                    f"Non-integer count at row {i}, species '{sp}': {value!r}"
                )
                continue
            if count < 0:
                errors.append(f"Negative count at row {i}, species '{sp}': {count}")

    only_tree, only_counts = check_tips_match_species(tree_file, species)
    if only_tree:
        errors.append(
            f"Species in tree but missing from count table: {sorted(only_tree)}"
        )
    if only_counts:
        errors.append(
            f"Species in count table but missing from tree: {sorted(only_counts)}"
        )

    return errors


def validate_output(output_dir: Path | str) -> list[str]:
    """Lightweight check that a CAFE-5 run produced results.

    CAFE-5 writes a ``<model>_results.txt`` (e.g. ``Base_results.txt`` /
    ``Gamma_results.txt``) plus per-family tables into its output directory.
    """
    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        return [f"CAFE-5 output directory not found: {output_dir}"]
    if not list(output_dir.glob("*_results.txt")):
        return [
            f"No CAFE-5 '*_results.txt' file found in {output_dir}; "
            "the run may have failed."
        ]
    return []


def build_command(
    count_file: Path | str,
    tree_file: Path | str,
    n_gamma_cats: int | None = None,
    extra_args: list[str] | str | None = None,
    tool_path: str = "cafe5",
) -> list[str]:
    """Build the CAFE-5 command: ``cafe5 -i <counts> -t <tree> [-k K] [extra]``."""
    cmd = [str(tool_path), "-i", str(count_file), "-t", str(tree_file)]
    if n_gamma_cats:
        cmd += ["-k", str(n_gamma_cats)]
    if extra_args:
        cmd += (
            shlex.split(extra_args)
            if isinstance(extra_args, str)
            else [str(arg) for arg in extra_args]
        )
    return cmd
