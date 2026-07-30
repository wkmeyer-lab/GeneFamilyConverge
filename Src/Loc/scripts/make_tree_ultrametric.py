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
- Dating uses a single penalized-likelihood fit at ``--smoothing`` (default
  100). ``--cross-validate`` auto-selects smoothing but scales ~O(taxa), so it
  is only feasible for small trees.
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
from convgeno.io.fasta import FASTA_EXTENSIONS
from convgeno.utils.command_runner import check_tool_available
from convgeno.utils.logging import setup_logging
from convgeno.validation.trees import (
    check_tips_match_species,
    is_ultrametric,
    ultrametric_deviation,
    validate_tree,
)

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
# User-supplied species tree helpers
# -------------------------------------------------------------------


def _species_from_dir(species_dir: Path) -> set[str]:
    """Return proteome species names (FASTA basenames) present in *species_dir*."""
    if not species_dir.is_dir():
        return set()
    return {
        entry.stem
        for entry in species_dir.iterdir()
        if entry.is_file() and entry.suffix.lower() in FASTA_EXTENSIONS
    }


def _report_tip_match(tree_path: Path, species_dir: Path | None) -> None:
    """Warn if the tree's tips don't match the proteome species set.

    This is the run-time re-validation of a user-supplied tree (the proteomes
    may not have existed when ``convgeno init`` recorded the tree). It only
    warns — the strict tip/count check lives in ``prepare_cafe_inputs.py``.
    """
    if species_dir is None:
        return
    species = _species_from_dir(species_dir)
    if not species:
        logger.info(
            "No proteome FASTAs in %s; deferring the tip-label check.", species_dir
        )
        return
    in_tree, in_proteomes = check_tips_match_species(tree_path, sorted(species))
    if in_proteomes:
        logger.warning(
            "%d proteome species are missing from the tree: %s",
            len(in_proteomes),
            sorted(in_proteomes)[:10],
        )
    if in_tree:
        logger.warning(
            "%d tree tips are not among the proteomes: %s",
            len(in_tree),
            sorted(in_tree)[:10],
        )
    if not in_tree and not in_proteomes:
        logger.info("Tree tips match all %d proteome species.", len(species))


def _run_assume_ultrametric(
    input_tree: Path | None,
    output: Path,
    species_dir: Path | None,
    stats_json: Path | None,
) -> int:
    """Validate a user tree as ultrametric and copy it to *output* (no r8s).

    Reads no OrthoFinder output. Returns a process exit code.
    """
    if input_tree is None:
        logger.error(
            "--assume-ultrametric requires --input-tree (or config species_tree.path)."
        )
        return 1
    if not input_tree.is_file():
        logger.error("Input species tree not found: %s", input_tree)
        return 1

    issues = validate_tree(input_tree)
    if issues:
        logger.error("Input tree is not usable: %s", "; ".join(issues))
        return 1
    if not is_ultrametric(input_tree):
        logger.error(
            "--assume-ultrametric was given, but %s is NOT ultrametric within "
            "tolerance (max-min root-to-tip deviation = %g). Re-run WITHOUT "
            "--assume-ultrametric so r8s can date it.",
            input_tree,
            ultrametric_deviation(input_tree),
        )
        return 1

    _report_tip_match(input_tree, species_dir)

    text = input_tree.read_text(encoding="utf-8")
    if not text.endswith("\n"):
        text += "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")

    print(f"ultrametric tree (user-supplied, r8s skipped) -> {output}")
    print("post-validation: tree is rooted, binary, and ultrametric.")

    if stats_json:
        stats_json.parent.mkdir(parents=True, exist_ok=True)
        with open(stats_json, "w", encoding="utf-8") as jf:
            json.dump(
                {
                    "mode": "assume_ultrametric",
                    "input_tree": str(input_tree),
                    "out_tree": str(output),
                },
                jf,
                indent=2,
            )
    return 0


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
        "--input-tree",
        type=Path,
        default=None,
        help="Date THIS rooted Newick instead of discovering OrthoFinder's "
        "SpeciesTree_rooted.txt (falls back to config species_tree.path).",
    )
    parser.add_argument(
        "--assume-ultrametric",
        action="store_true",
        help="Treat --input-tree as already ultrametric: validate it and copy it "
        "to the output path WITHOUT running r8s (falls back to config "
        "species_tree.is_ultrametric).",
    )
    parser.add_argument(
        "--species-dir",
        type=Path,
        default=None,
        help="Directory of proteome FASTAs; when given, the tree's tips are "
        "checked against these species (falls back to config "
        "orthofinder.input_dir).",
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
        "--smoothing",
        type=float,
        default=None,
        help="PL smoothing for a single fit "
        "(default: config ultrametric.smoothing, else 100).",
    )
    parser.add_argument(
        "--cross-validate",
        action="store_true",
        help="Cross-validate smoothing instead of a single fit. Accurate but "
        "scales ~O(taxa) — only practical for small trees.",
    )
    parser.add_argument(
        "--root-age",
        type=float,
        default=None,
        help="Root age used when no calibration is given (relative-time tree; "
        "default: config ultrametric.root_age, else 1).",
    )
    parser.add_argument(
        "--skip-if-unavailable",
        action="store_true",
        help="If r8s is not installed, log a message and exit 0 instead of "
        "failing (for automatic post-OrthoFinder runs).",
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

    output = args.output or _dig(config, "ultrametric", "output")
    if output is None:
        raise ValueError(
            "Output path not given: pass -o/--output or set ultrametric.output "
            "in --config."
        )

    input_tree = args.input_tree or _dig(config, "species_tree", "path")
    assume_ultrametric = args.assume_ultrametric or bool(
        _dig(config, "species_tree", "is_ultrametric")
    )
    species_dir = args.species_dir or _dig(config, "orthofinder", "input_dir")

    nsites = args.nsites
    if nsites is None:
        cfg_nsites = _dig(config, "ultrametric", "nsites")
        if isinstance(cfg_nsites, int):
            nsites = cfg_nsites  # a str like "auto" means auto-derive
    if nsites is None:
        cfg_st_nsites = _dig(config, "species_tree", "num_sites")
        if isinstance(cfg_st_nsites, int):
            nsites = cfg_st_nsites

    calibrations = [_parse_calibration(spec) for spec in args.calibration]
    calibrations.extend(_pairs_to_calibrations(args.pairs, args.cal_ages))
    if not calibrations:
        cfg_cals = _dig(config, "ultrametric", "calibrations")
        if cfg_cals:
            calibrations = _calibrations_from_config(cfg_cals)
    # No calibration is allowed: the tree is anchored at root_age (relative time).

    r8s_path = args.r8s_path or _dig(tools, "r8s", "command") or "r8s"

    root_age = args.root_age
    if root_age is None:
        cfg_root_age = _dig(config, "ultrametric", "root_age")
        root_age = float(cfg_root_age) if cfg_root_age is not None else 1.0

    smoothing = args.smoothing
    if smoothing is None:
        cfg_smoothing = _dig(config, "ultrametric", "smoothing")
        smoothing = float(cfg_smoothing) if cfg_smoothing is not None else 100.0

    cross_validate = args.cross_validate or bool(
        _dig(config, "ultrametric", "cross_validate")
    )

    # OrthoFinder output is needed unless the tree is taken as-is
    # (--assume-ultrametric) or a user tree + explicit nsites make it redundant.
    needs_of_output = not assume_ultrametric and not (
        input_tree is not None and nsites is not None
    )
    if needs_of_output and of_dir is None:
        raise ValueError(
            "OrthoFinder output dir not given: pass it positionally or set "
            "orthofinder.output_dir in --config. It is needed to locate the "
            "species tree and/or the alignment used for nsites."
        )

    return {
        "of_dir": Path(of_dir) if of_dir is not None else None,
        "output": Path(output),
        "input_tree": Path(input_tree) if input_tree is not None else None,
        "assume_ultrametric": assume_ultrametric,
        "species_dir": Path(species_dir) if species_dir is not None else None,
        "nsites": nsites,
        "calibrations": calibrations,
        "r8s_path": r8s_path,
        "smoothing": smoothing,
        "cross_validate": cross_validate,
        "root_age": root_age,
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
        input_tree = resolved["input_tree"]
        work_dir = args.work_dir or (output.parent / "r8s_work")

        # A user-supplied tree that is already ultrametric skips r8s entirely
        # (and reads no OrthoFinder output).
        if resolved["assume_ultrametric"]:
            return _run_assume_ultrametric(
                input_tree, output, resolved["species_dir"], args.stats_json
            )

        if (
            args.skip_if_unavailable
            and not args.dry_run
            and not check_tool_available(resolved["r8s_path"])
        ):
            logger.warning(
                "r8s ('%s') is not installed; skipping the ultrametric step "
                "(--skip-if-unavailable). Build it with tools/r8s/install_r8s.sh.",
                resolved["r8s_path"],
            )
            print("r8s not installed; ultrametric step skipped.")
            return 0

        # Re-validate a user tree's tips against the proteomes at run time.
        if input_tree is not None:
            _report_tip_match(input_tree, resolved["species_dir"])

        stats = make_ultrametric(
            resolved["of_dir"],
            resolved["calibrations"],
            out_tree=output,
            work_dir=work_dir,
            input_tree=input_tree,
            nsites=resolved["nsites"],
            r8s_path=resolved["r8s_path"],
            method=args.method,
            algorithm=args.algorithm,
            smoothing=resolved["smoothing"],
            cross_validate=resolved["cross_validate"],
            root_age=resolved["root_age"],
            dry_run=args.dry_run,
        )
    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        logger.error("%s", exc)
        return 1

    dating = (
        "cross-validated smoothing"
        if stats.get("cross_validate")
        else f"fixed smoothing={stats['smoothing']:g}"
    )
    calib = (
        "none (relative-time, root-anchored)"
        if stats.get("relative_time")
        else str(stats["n_calibrations"])
    )
    print(f"nsites:        {stats['nsites']}")
    print(f"calibrations:  {calib}")
    print(f"tips ({len(stats['tips'])}):     {', '.join(stats['tips'])}")
    print(f"dating:        {args.method} ({dating})")
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
