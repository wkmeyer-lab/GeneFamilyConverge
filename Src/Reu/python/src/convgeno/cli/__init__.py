"""Entry point for the convgeno command-line interface."""

from __future__ import annotations

import argparse
import sys


def main() -> None:
    """Entry point for the convgeno command-line interface."""
    parser = argparse.ArgumentParser(
        prog="convgeno",
        description=(
            "convgeno: a pipeline for gene family copy number analysis "
            "across species."
        ),
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = False

    init_parser = subparsers.add_parser(
        "init",
        help=(
            "Interactive first-time setup. Detects SLURM partitions and "
            "writes pipeline_config.yaml."
        ),
    )
    init_parser.add_argument(
        "--output",
        default="pipeline_config.yaml",
        help="Path to write the config file (default: pipeline_config.yaml)",
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "init":
        from convgeno.cli.init_cmd import run_init

        run_init(output_path=args.output)
