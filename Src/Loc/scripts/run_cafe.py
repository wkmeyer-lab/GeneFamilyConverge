#!/usr/bin/env python3
"""
run_cafe.py — Run CAFE-5 to model gene-family expansion/contraction.

Given the CAFE-5 count matrix and the rooted/binary/ultrametric species tree
(both produced by ``prepare_cafe_inputs.py``), this runs CAFE-5, optionally with
a categorical phenotype tree (the ``-y`` multi-lambda tree) so gain/loss rates
can differ between trait states.

This is the final step of the ``convgeno run`` DAG (``workflow/Snakefile``); it is
invoked with explicit paths::

    python Src/Loc/scripts/run_cafe.py \\
        --counts Data/interim/cafe_input/cafe_input.tsv \\
        --tree   Data/interim/cafe_input/species_tree.nwk \\
        [--lambda-tree Data/processed/orthofinder_.../lambda_tree.nwk] \\
        --output-dir Data/processed/cafe_results

CAFE-5 is NOT on conda; install it separately and put it on PATH (or pass
``--cafe-path``). See the README's CAFE-5 section.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from convgeno.external import cafe
from convgeno.utils.command_runner import CommandError, check_tool_available
from convgeno.utils.command_runner import run as run_command
from convgeno.utils.logging import setup_logging

logger = logging.getLogger(__name__)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_cafe",
        description="Run CAFE-5 on a gene-family count matrix + ultrametric tree.",
    )
    parser.add_argument(
        "--counts",
        type=Path,
        required=True,
        help="CAFE-5 count matrix (cafe_input.tsv from prepare_cafe_inputs.py).",
    )
    parser.add_argument(
        "--tree",
        type=Path,
        required=True,
        help="Rooted, binary, ultrametric species tree (Newick).",
    )
    parser.add_argument(
        "--lambda-tree",
        type=Path,
        default=None,
        help=(
            "Optional categorical phenotype tree (CAFE -y). A missing file is "
            "skipped with a warning (single-lambda run); a malformed one is an error."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="CAFE-5 output directory (-o).",
    )
    parser.add_argument(
        "-k",
        "--gamma-categories",
        type=int,
        default=None,
        help="Number of gamma rate categories (CAFE -k). Omit for the base model.",
    )
    parser.add_argument(
        "-p",
        "--poisson",
        action="store_true",
        help="Use a Poisson root-frequency distribution (CAFE -p).",
    )
    parser.add_argument(
        "--cafe-path",
        default="cafe5",
        help="Path/name of the CAFE-5 executable (default: cafe5).",
    )
    parser.add_argument(
        "--extra",
        default=None,
        help="Extra arguments passed verbatim to cafe5 (quoted string).",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )
    parser.add_argument("--log-file", type=Path, default=None)
    return parser


def _resolve_lambda_tree(lambda_tree: Path | None, tree: Path) -> Path | None:
    """Validate the optional -y tree. Return it, or ``None`` to run single-lambda.

    A missing file is treated as "phenotype step skipped" (warn, drop -y). A file
    that exists but fails validation raises ``ValueError`` (a real problem).
    """
    if lambda_tree is None:
        return None
    if not lambda_tree.is_file():
        logger.warning(
            "Phenotype (-y) tree %s not found; running CAFE-5 with a single "
            "lambda. (The categorical phenotype tree step may have been skipped.)",
            lambda_tree,
        )
        return None
    errors = cafe.validate_lambda_tree(lambda_tree, tree)
    if errors:
        raise ValueError(
            "Categorical phenotype (-y) tree failed validation:\n  - "
            + "\n  - ".join(errors)
        )
    return lambda_tree


def run_cafe(
    counts: Path,
    tree: Path,
    output_dir: Path,
    *,
    lambda_tree: Path | None = None,
    gamma_categories: int | None = None,
    poisson: bool = False,
    cafe_path: str = "cafe5",
    extra: str | None = None,
) -> int:
    """Validate inputs, run CAFE-5, verify outputs. Returns a process exit code."""
    if not check_tool_available(cafe_path):
        logger.error(
            "CAFE-5 executable '%s' not found on PATH. CAFE-5 is not on conda; "
            "install it separately and put it on PATH (or pass --cafe-path). "
            "See the README's CAFE-5 section.",
            cafe_path,
        )
        return 1

    input_errors = cafe.validate_input(counts, tree)
    if input_errors:
        logger.error("CAFE-5 input validation failed:")
        for msg in input_errors:
            logger.error("  - %s", msg)
        return 1

    try:
        effective_lambda = _resolve_lambda_tree(lambda_tree, tree)
    except ValueError as exc:
        logger.error("%s", exc)
        return 1

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # build_command does not set -o; CAFE-5 otherwise writes to ./results.
    extra_args: list[str] = ["-o", str(output_dir)]
    if poisson:
        extra_args.append("-p")
    if extra:
        import shlex

        extra_args.extend(shlex.split(extra))

    cmd = cafe.build_command(
        counts,
        tree,
        n_gamma_cats=gamma_categories,
        extra_args=extra_args,
        tool_path=cafe_path,
        lambda_tree=effective_lambda,
    )

    try:
        run_command(cmd, log_dir=output_dir)
    except CommandError as exc:
        logger.error("CAFE-5 failed: %s", exc)
        return 1

    output_errors = cafe.validate_output(output_dir)
    if output_errors:
        for msg in output_errors:
            logger.error("  - %s", msg)
        return 1

    results = sorted(p.name for p in output_dir.glob("*_results.txt"))
    print(f"CAFE-5 complete. Results in {output_dir}: {', '.join(results)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run CAFE-5. Returns an exit code."""
    args = _build_parser().parse_args(argv)
    setup_logging(level=args.log_level, log_file=args.log_file)
    return run_cafe(
        counts=args.counts,
        tree=args.tree,
        output_dir=args.output_dir,
        lambda_tree=args.lambda_tree,
        gamma_categories=args.gamma_categories,
        poisson=args.poisson,
        cafe_path=args.cafe_path,
        extra=args.extra,
    )


if __name__ == "__main__":
    sys.exit(main())
