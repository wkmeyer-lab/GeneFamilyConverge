"""
convgeno.validation.trees
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Validate phylogenetic tree inputs for CAFE-5, which requires a **rooted,
strictly binary, ultrametric** Newick tree whose tip labels match the species
in the gene-count table.

Uses ``Bio.Phylo`` (biopython is already a dependency) — no extra tree library.
"""

from __future__ import annotations

import logging
from io import StringIO
from pathlib import Path

from Bio import Phylo

logger = logging.getLogger(__name__)

__all__ = [
    "validate_tree",
    "check_tips_match_species",
    "is_ultrametric",
    "root_to_tip_depths",
    "ultrametric_deviation",
]


def _load_tree(source: Path | str):
    """Load a tree from a path or a raw Newick string."""
    if isinstance(source, Path):
        text = source.read_text(encoding="utf-8")
    else:
        text = str(source)
        if "(" not in text and ";" not in text and Path(text).is_file():
            text = Path(text).read_text(encoding="utf-8")
    return Phylo.read(StringIO(text), "newick")


def _tip_depths(tree) -> dict:
    """Map each terminal clade to its root-to-tip distance (branch-length sum)."""
    depths: dict = {}

    def _recurse(clade, acc: float) -> None:
        acc += clade.branch_length or 0.0
        if clade.is_terminal():
            depths[clade] = acc
        for child in clade.clades:
            _recurse(child, acc)

    _recurse(tree.root, 0.0)
    return depths


def validate_tree(path: Path | str) -> list[str]:
    """Check a Newick tree is parseable, rooted, and strictly binary.

    Parameters
    ----------
    path : Path or str
        A path to a Newick file, or the Newick string itself.

    Returns
    -------
    list[str]
        Error messages; an empty list means the tree is valid.
    """
    try:
        tree = _load_tree(path)
    except Exception as exc:  # noqa: BLE001 - report any parse failure uniformly
        return [f"Could not parse tree: {exc}"]

    errors: list[str] = []

    root = tree.root
    if len(root.clades) != 2:
        errors.append(
            f"Tree is not rooted as a bifurcating tree: the root has "
            f"{len(root.clades)} children (expected 2)."
        )

    polytomies = sum(1 for nt in tree.get_nonterminals() if len(nt.clades) > 2)
    if polytomies:
        errors.append(
            f"Tree contains {polytomies} polytomy/polytomies (internal node with "
            f">2 children); CAFE-5 requires a strictly binary tree."
        )

    return errors


def check_tips_match_species(
    tree_path: Path | str, species_list: list[str]
) -> tuple[set, set]:
    """Compare tree tip labels against a species list.

    Returns
    -------
    tuple[set, set]
        ``(in_tree_not_in_list, in_list_not_in_tree)``.  Both empty means the
        tip labels and the species list agree exactly.
    """
    tree = _load_tree(tree_path)
    tips = {tip.name for tip in tree.get_terminals()}
    species = set(species_list)
    return (tips - species, species - tips)


def root_to_tip_depths(tree_path: Path | str) -> dict[str, float]:
    """Map each tip label to its root-to-tip distance (sum of branch lengths).

    Reuses the same ``Bio.Phylo`` loading as the other validators. Useful for
    *reporting* why a tree is or isn't ultrametric (equal depths == ultrametric).
    """
    tree = _load_tree(tree_path)
    return {clade.name: depth for clade, depth in _tip_depths(tree).items()}


def ultrametric_deviation(tree_path: Path | str) -> float:
    """Return ``max_depth - min_depth`` across all root-to-tip paths.

    ``0.0`` for a perfectly ultrametric tree; larger values mean the tips are
    increasingly out of alignment. Returns ``0.0`` for a tree with no tips.
    """
    depths = list(_tip_depths(_load_tree(tree_path)).values())
    if not depths:
        return 0.0
    return max(depths) - min(depths)


def is_ultrametric(tree_path: Path | str, tolerance: float = 1e-3) -> bool:
    """Return True if every root-to-tip path has (nearly) the same length.

    ``tolerance`` is interpreted **relative to the tree height**: the tree is
    ultrametric if ``max_depth - min_depth <= tolerance * max_depth``.  A
    relative tolerance is scale-invariant, so it works whether branch lengths
    are absolute times (Myr) or relative units.
    """
    tree = _load_tree(tree_path)
    depths = _tip_depths(tree)
    if not depths:
        return False
    values = list(depths.values())
    lo, hi = min(values), max(values)
    scale = hi if hi > 0 else 1.0
    return (hi - lo) <= tolerance * scale
