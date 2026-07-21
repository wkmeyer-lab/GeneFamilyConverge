#!/usr/bin/env python3
"""
Make an OrthoFinder species tree ultrametric with r8s (for CAFE-5 input).

CAFE-5 needs a rooted, ultrametric time tree, but OrthoFinder emits a
non-ultrametric ``SpeciesTree_rooted.txt``.  This script derives ``nsites`` from
the concatenated species-tree alignment, runs r8s with one or more fossil
calibrations, and writes the dated (ultrametric) Newick.

Basic usage (point calibration, tutorial-style)::

    python Src/Loc/scripts/make_tree_ultrametric.py \\
        Data/processed/orthofinder_single_<ts> \\
        -o Data/interim/cafe_input/species_tree_ultrametric.nwk \\
        -p 'human,cat' -c 94

Named calibrations (repeatable), point or age window::

    python Src/Loc/scripts/make_tree_ultrametric.py <of_output_dir> -o tree.nwk \\
        --calibration humancat:human,cat:94 \\
        --calibration root:human,dog:100-120

Notes
-----
- Calibration taxa must be tip labels in ``SpeciesTree_rooted.txt`` (i.e. the
  OrthoFinder proteome basenames).
- ``nsites`` is auto-derived from ``MultipleSequenceAlignments/
  SpeciesTreeAlignment.fa`` (MSA mode); override with ``--nsites`` if needed.
- r8s is not on conda; install it separately and pass ``--r8s-path`` (or put it
  on ``PATH``).  Use ``--dry-run`` to write the control file without running r8s.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path

from convgeno.external.r8s import Calibration, make_ultrametric
from convgeno.utils.logging import setup_logging

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------
# Calibration parsing
# -------------------------------------------------------------------


def _sanitize_node_name(name: str) -> str:
    """Reduce an MRCA node label to r8s-safe characters."""
    cleaned = re.sub(r"[^A-Za-z0-9_]", "_", name.strip())
    return cleaned or "node"


def _is_number(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


def _age_kwargs(age_spec: str) -> dict:
    """Parse an age token into Calibration kwargs.

    ``"94"`` -> fixed age; ``"100-120"`` -> window; ``"100-"`` / ``"-120"`` ->
    one-sided window.
    """
    token = age_spec.strip()
    if _is_number(token):
        return {"age": float(token)}
    if "-" in token:
        low, _, high = token.partition("-")
        min_age = float(low) if low.strip() else None
        max_age = float(high) if high.strip() else None
        if min_age is None and max_age is None:
            raise ValueError(f"Invalid age window: {age_spec!r}")
        return {"min_age": min_age, "max_age": max_age}
    raise ValueError(f"Could not parse age spec: {age_spec!r} (use AGE or MIN-MAX)")


def _parse_calibration(spec: str) -> Calibration:
    """Parse ``NAME:SP1,SP2:AGESPEC`` into a Calibration."""
    parts = spec.split(":")
    if len(parts) != 3:
        raise ValueError(f"--calibration must be 'NAME:SP1,SP2:AGE' (got {spec!r})")
    name, taxa_str, age_spec = parts
    taxa = tuple(t.strip() for t in taxa_str.split(","))
    if len(taxa) != 2:
        raise ValueError(f"Calibration taxa must be 'SP1,SP2' (got {taxa_str!r})")
    return Calibration(
        name=_sanitize_node_name(name), taxa=taxa, **_age_kwargs(age_spec)
    )


def _pairs_to_calibrations(pairs: list[str], ages: list[str]) -> list[Calibration]:
    """Convert tutorial-style ``-p``/``-c`` arguments into Calibrations."""
    if len(pairs) != len(ages):
        raise ValueError(
            f"Each -p/--pair needs a matching -c/--cal: got {len(pairs)} pair(s) "
            f"and {len(ages)} age(s)."
        )
    cals: list[Calibration] = []
    for pair, age in zip(pairs, ages, strict=True):
        taxa = tuple(t.strip() for t in pair.split(","))
        if len(taxa) != 2:
            raise ValueError(f"-p/--pair must be 'SP1,SP2' (got {pair!r})")
        name = _sanitize_node_name(f"{taxa[0]}_{taxa[1]}")
        cals.append(Calibration(name=name, taxa=taxa, **_age_kwargs(age)))
    return cals


def _collect_calibrations(args: argparse.Namespace) -> list[Calibration]:
    cals = [_parse_calibration(spec) for spec in args.calibration]
    cals.extend(_pairs_to_calibrations(args.pairs, args.cal_ages))
    if not cals:
        raise ValueError(
            "Provide at least one calibration via --calibration "
            "or -p/--pair + -c/--cal."
        )
    return cals


# -------------------------------------------------------------------
# Argument parsing
# -------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="make_tree_ultrametric",
        description="Make an OrthoFinder species tree ultrametric with r8s.",
    )
    parser.add_argument(
        "orthofinder_output_dir",
        type=Path,
        help="OrthoFinder -o output dir (or a Results_* dir directly).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output path for the ultrametric Newick tree.",
    )
    parser.add_argument(
        "--calibration",
        action="append",
        default=[],
        metavar="NAME:SP1,SP2:AGE",
        help="Named calibration; AGE is a number or MIN-MAX window. Repeatable.",
    )
    parser.add_argument(
        "-p",
        "--pair",
        action="append",
        default=[],
        dest="pairs",
        metavar="'SP1,SP2'",
        help="Tutorial-style calibrated taxon pair (repeatable, paired with -c).",
    )
    parser.add_argument(
        "-c",
        "--cal",
        action="append",
        default=[],
        dest="cal_ages",
        metavar="AGE",
        help="Age for the corresponding -p pair (repeatable).",
    )
    parser.add_argument(
        "--nsites",
        type=int,
        default=None,
        help="Override the alignment column count (default: auto from the MSA).",
    )
    parser.add_argument(
        "--r8s-path",
        default="r8s",
        help="Path to the r8s binary (default: 'r8s' on PATH).",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Scratch dir for the control file / r8s output "
        "(default: <output-dir>/r8s_work).",
    )
    parser.add_argument(
        "--method",
        default="pl",
        help="r8s divtime method (default: pl).",
    )
    parser.add_argument(
        "--algorithm",
        default="tn",
        help="r8s divtime algorithm (default: tn).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write the r8s control file and log the command without running r8s.",
    )
    parser.add_argument(
        "--stats-json",
        type=Path,
        default=None,
        help="Write the run's summary stats to this JSON file.",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Also write log messages to this file.",
    )
    return parser


# -------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run the r8s ultrametric step.  Returns exit code."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    setup_logging(level=args.log_level, log_file=args.log_file)

    work_dir = args.work_dir or (args.output.parent / "r8s_work")

    try:
        calibrations = _collect_calibrations(args)
        stats = make_ultrametric(
            args.orthofinder_output_dir,
            calibrations,
            out_tree=args.output,
            work_dir=work_dir,
            nsites=args.nsites,
            r8s_path=args.r8s_path,
            method=args.method,
            algorithm=args.algorithm,
            dry_run=args.dry_run,
        )
    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        logger.error("%s", exc)
        return 1

    print(f"nsites:        {stats['nsites']}")
    print(f"calibrations:  {stats['n_calibrations']}")
    print(f"tips ({len(stats['tips'])}):     {', '.join(stats['tips'])}")
    print(f"control file:  {stats['r8s_ctl']}")
    if stats.get("dry_run"):
        print("dry run: r8s was not executed; no tree written.")
    else:
        print(f"ultrametric tree -> {stats['out_tree']}")

    if args.stats_json:
        args.stats_json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.stats_json, "w", encoding="utf-8") as jf:
            json.dump(stats, jf, indent=2)
        logger.info("Wrote stats to %s", args.stats_json)

    return 0


if __name__ == "__main__":
    sys.exit(main())
