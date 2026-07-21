#!/usr/bin/env python3
"""
Make an OrthoFinder species tree ultrametric with r8s (for CAFE-5 input).

CAFE-5 needs a rooted, ultrametric time tree, but OrthoFinder emits a
non-ultrametric ``SpeciesTree_rooted.txt``.  This script derives ``nsites`` from
the concatenated species-tree alignment, runs r8s with one or more fossil
calibrations, writes the dated (ultrametric) Newick, and post-validates it.

Config-driven (used by the workflows)::

    python Src/Loc/scripts/make_tree_ultrametric.py \\
        --config Src/Loc/configs/example_config.yaml \\
        --tools  Src/Loc/configs/tool_paths.yaml

Explicit, point calibration (tutorial-style)::

    python Src/Loc/scripts/make_tree_ultrametric.py <of_output_dir> \\
        -o Data/interim/cafe_input/species_tree_ultrametric.nwk \\
        -p 'human,cat' -c 94

Named calibrations (repeatable), point or age window::

    python Src/Loc/scripts/make_tree_ultrametric.py <of_output_dir> -o tree.nwk \\
        --calibration humancat:human,cat:94 \\
        --calibration root:human,dog:100-120

Notes
-----
- Calibration taxa must be tip labels in ``SpeciesTree_rooted.txt`` (the
  OrthoFinder proteome basenames).
- ``nsites`` is auto-derived from ``MultipleSequenceAlignments/
  SpeciesTreeAlignment.fa`` (MSA mode); override with ``--nsites``.
- r8s is not on conda; install it separately and set ``r8s.command`` in
  ``tool_paths.yaml`` or pass ``--r8s-path``.  ``--dry-run`` writes the control
  file without running r8s.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path

import yaml

from convgeno.external.r8s import Calibration, make_ultrametric
from convgeno.utils.logging import setup_logging
from convgeno.validation.trees import is_ultrametric, validate_tree

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


def _calibrations_from_config(entries: list[dict]) -> list[Calibration]:
    """Build Calibrations from the config's ``ultrametric.calibrations`` list."""
    cals: list[Calibration] = []
    for entry in entries:
        taxa = tuple(entry["taxa"])
        if len(taxa) != 2:
            raise ValueError(f"Config calibration taxa must be a pair (got {taxa!r})")
        name = entry.get("name") or f"{taxa[0]}_{taxa[1]}"
        kwargs = {
            k: entry[k]
            for k in ("age", "min_age", "max_age")
            if entry.get(k) is not None
        }
        if not kwargs:
            raise ValueError(
                f"Config calibration {name!r} needs 'age' or 'min_age'/'max_age'."
            )
        cals.append(Calibration(name=_sanitize_node_name(name), taxa=taxa, **kwargs))
    return cals


# -------------------------------------------------------------------
# Config helpers
# -------------------------------------------------------------------


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
        nargs="?",
        default=None,
        help="OrthoFinder -o output dir (or a Results_* dir). "
        "Falls back to config orthofinder.output_dir.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output Newick path. Falls back to config ultrametric.output.",
    )
    parser.add_argument(
        "--config", type=Path, default=None, help="Pipeline config YAML."
    )
    parser.add_argument("--tools", type=Path, default=None, help="tool_paths.yaml.")
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
        help="Override the alignment column count (default: auto / config).",
    )
    parser.add_argument(
        "--r8s-path",
        default=None,
        help="Path to r8s (default: tool_paths r8s.command, else 'r8s').",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Scratch dir (default: <output-dir>/r8s_work).",
    )
    parser.add_argument(
        "--method", default="pl", help="r8s divtime method (default: pl)."
    )
    parser.add_argument(
        "--algorithm", default="tn", help="r8s divtime algorithm (default: tn)."
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


def _resolve_inputs(args: argparse.Namespace) -> dict:
    """Merge CLI args with config/tools files into resolved run parameters."""
    config = _load_yaml(args.config)
    tools = _load_yaml(args.tools)

    of_dir = args.orthofinder_output_dir or _dig(config, "orthofinder", "output_dir")
    if of_dir is None:
        raise ValueError(
            "OrthoFinder output dir not given: pass it positionally or set "
            "orthofinder.output_dir in --config."
        )

    output = args.output or _dig(config, "ultrametric", "output")
    if output is None:
        raise ValueError(
            "Output path not given: pass -o/--output or set ultrametric.output "
            "in --config."
        )

    nsites = args.nsites
    if nsites is None:
        cfg_nsites = _dig(config, "ultrametric", "nsites")
        if isinstance(cfg_nsites, int):
            nsites = cfg_nsites  # a str like "auto" means auto-derive

    calibrations = [_parse_calibration(spec) for spec in args.calibration]
    calibrations.extend(_pairs_to_calibrations(args.pairs, args.cal_ages))
    if not calibrations:
        cfg_cals = _dig(config, "ultrametric", "calibrations")
        if cfg_cals:
            calibrations = _calibrations_from_config(cfg_cals)
    if not calibrations:
        raise ValueError(
            "No calibrations: use --calibration / -p+-c, or add "
            "ultrametric.calibrations to --config."
        )

    r8s_path = args.r8s_path or _dig(tools, "r8s", "command") or "r8s"

    return {
        "of_dir": Path(of_dir),
        "output": Path(output),
        "nsites": nsites,
        "calibrations": calibrations,
        "r8s_path": r8s_path,
    }


# -------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run the r8s ultrametric step.  Returns exit code."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    setup_logging(level=args.log_level, log_file=args.log_file)

    try:
        resolved = _resolve_inputs(args)
        output = resolved["output"]
        work_dir = args.work_dir or (output.parent / "r8s_work")

        stats = make_ultrametric(
            resolved["of_dir"],
            resolved["calibrations"],
            out_tree=output,
            work_dir=work_dir,
            nsites=resolved["nsites"],
            r8s_path=resolved["r8s_path"],
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
        issues = validate_tree(output)
        if not is_ultrametric(output):
            issues.append("tree is not ultrametric within tolerance")
        if issues:
            logger.warning(
                "Post-validation found %d issue(s) in %s", len(issues), output
            )
            print("post-validation WARNINGS:")
            for msg in issues:
                print(f"  ! {msg}")
        else:
            print("post-validation: tree is rooted, binary, and ultrametric.")

    if args.stats_json:
        args.stats_json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.stats_json, "w", encoding="utf-8") as jf:
            json.dump(stats, jf, indent=2)
        logger.info("Wrote stats to %s", args.stats_json)

    return 0


if __name__ == "__main__":
    sys.exit(main())
