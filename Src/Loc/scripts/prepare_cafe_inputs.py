#!/usr/bin/env python3
"""
Convert OrthoFinder output into CAFE-5 inputs.

CAFE-5 requires:
    1. A tab-delimited count file: ``Desc<TAB>Family ID<TAB><species...>``.
    2. A rooted, binary, ultrametric species tree (from make_tree_ultrametric).

This script:
    - locates ``Orthogroups_GeneCount.tsv`` in the OrthoFinder output,
    - reformats it into CAFE-5's layout (Desc + Family ID + species counts),
    - validates the count file against the ultrametric tree (tips match,
      counts are non-negative integers, tree is binary/rooted/ultrametric),
    - copies both into the CAFE-5 input directory.

Config-driven (used by the workflows)::

    python Src/Loc/scripts/prepare_cafe_inputs.py \\
        --config Src/Loc/configs/example_config.yaml

Explicit paths::

    python Src/Loc/scripts/prepare_cafe_inputs.py \\
        --orthofinder-dir Data/processed/orthofinder_single_<ts> \\
        --tree Data/interim/cafe_input/species_tree_ultrametric.nwk \\
        --output-dir Data/interim/cafe_input
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path

import yaml

from convgeno.external.cafe import format_gene_counts_for_cafe, validate_input
from convgeno.external.r8s import find_results_dir
from convgeno.utils.logging import setup_logging

logger = logging.getLogger(__name__)

GENE_COUNT_RELPATH = Path("Orthogroups") / "Orthogroups_GeneCount.tsv"


def _load_yaml(path: Path | None) -> dict:
    if path is None:
        return {}
    return yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}


def _dig(mapping: dict, *keys: str):
    node = mapping
    for key in keys:
        if not isinstance(node, dict):
            return None
        node = node.get(key)
    return node


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="prepare_cafe_inputs",
        description="Reformat OrthoFinder gene counts + validate tree for CAFE-5.",
    )
    parser.add_argument(
        "--config", type=Path, default=None, help="Pipeline config YAML."
    )
    parser.add_argument(
        "--orthofinder-dir",
        type=Path,
        default=None,
        help="OrthoFinder -o output dir (default: config orthofinder.output_dir).",
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=None,
        help="Explicit Orthogroups_GeneCount.tsv (default: found under the OF dir).",
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=None,
        help="Ultrametric species tree (default: config ultrametric.output).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="CAFE-5 input dir (default: <config interim_dir>/cafe_input "
        "or Data/interim/cafe_input).",
    )
    parser.add_argument(
        "--description",
        default="(null)",
        help="Value for the CAFE-5 'Desc' column (default: (null)).",
    )
    parser.add_argument(
        "--max-family-size",
        type=int,
        default=None,
        help="Exclude families with >= this count in any species (e.g. 100).",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )
    parser.add_argument("--log-file", type=Path, default=None)
    return parser


def _resolve(args: argparse.Namespace) -> dict:
    config = _load_yaml(args.config)

    of_dir = args.orthofinder_dir or _dig(config, "orthofinder", "output_dir")

    counts = args.counts
    if counts is None:
        if of_dir is None:
            raise ValueError(
                "Provide --counts, or --orthofinder-dir / config "
                "orthofinder.output_dir to locate Orthogroups_GeneCount.tsv."
            )
        counts = find_results_dir(of_dir) / GENE_COUNT_RELPATH

    tree = args.tree or _dig(config, "ultrametric", "output")
    if tree is None:
        raise ValueError(
            "Provide --tree, or set ultrametric.output in --config (the "
            "ultrametric tree from make_tree_ultrametric)."
        )

    if args.output_dir is not None:
        output_dir = args.output_dir
    else:
        interim = _dig(config, "outputs", "interim_dir") or "Data/interim"
        output_dir = Path(interim) / "cafe_input"

    # CLI overrides config; config's cafe.max_family_size is the default.
    max_family_size = args.max_family_size
    if max_family_size is None:
        max_family_size = _dig(config, "cafe", "max_family_size")

    return {
        "counts": Path(counts),
        "tree": Path(tree),
        "output_dir": Path(output_dir),
        "max_family_size": max_family_size,
    }


def main(argv: list[str] | None = None) -> int:
    """Prepare CAFE-5 inputs.  Returns exit code."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    setup_logging(level=args.log_level, log_file=args.log_file)

    try:
        resolved = _resolve(args)
        output_dir = resolved["output_dir"]
        count_out = output_dir / "cafe_input.tsv"
        tree_out = output_dir / "species_tree.nwk"

        stats = format_gene_counts_for_cafe(
            resolved["counts"],
            count_out,
            description=args.description,
            max_family_size=resolved["max_family_size"],
        )

        errors = validate_input(count_out, resolved["tree"])
        if errors:
            logger.error("CAFE-5 input validation failed:")
            for msg in errors:
                logger.error("  - %s", msg)
            return 1

        tree_out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(resolved["tree"], tree_out)
    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        logger.error("%s", exc)
        return 1

    print(f"species ({stats['n_species']}): {', '.join(stats['species'])}")
    print(f"families written: {stats['n_written']} / {stats['n_families']}")
    if "n_large_excluded" in stats:
        print(
            f"large families excluded: {stats['n_large_excluded']} "
            f"-> {stats['large_family_file']}"
        )
    print(f"count file -> {count_out}")
    print(f"species tree -> {tree_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
