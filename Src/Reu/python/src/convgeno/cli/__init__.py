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

    of_parser = subparsers.add_parser(
        "orthofinder",
        help="Generate and submit OrthoFinder SLURM jobs.",
    )
    of_subparsers = of_parser.add_subparsers(dest="orthofinder_command")

    of_gen = of_subparsers.add_parser(
        "generate",
        help="Generate the OrthoFinder SLURM script without submitting.",
    )
    of_gen.add_argument(
        "--config",
        default="pipeline_config.yaml",
        help="Path to pipeline config file (default: pipeline_config.yaml)",
    )
    of_gen.add_argument(
        "--script",
        default="slurm_scripts/orthofinder.sh",
        help=(
            "Path to write the generated SLURM script "
            "(default: slurm_scripts/orthofinder.sh)"
        ),
    )
    of_gen.add_argument(
        "--multinode",
        action="store_true",
        help=(
            "Generate three scripts for multi-node execution "
            "(prepare, search array, resume) instead of one."
        ),
    )
    of_gen.add_argument(
        "--script-dir",
        default="slurm_scripts",
        help=(
            "Directory to write generated scripts (used with --multinode, "
            "default: slurm_scripts)"
        ),
    )

    of_run = of_subparsers.add_parser(
        "run",
        help="Generate the OrthoFinder SLURM script and submit it.",
    )
    of_run.add_argument(
        "--config",
        default="pipeline_config.yaml",
        help="Path to pipeline config file (default: pipeline_config.yaml)",
    )
    of_run.add_argument(
        "--script",
        default="slurm_scripts/orthofinder.sh",
        help=(
            "Path to write the generated SLURM script "
            "(default: slurm_scripts/orthofinder.sh)"
        ),
    )
    of_run.add_argument(
        "-y", "--yes",
        action="store_true",
        help="Skip confirmation prompt and submit immediately.",
    )
    of_run.add_argument(
        "--multinode",
        action="store_true",
        help=(
            "Generate and submit three scripts for multi-node execution "
            "(prepare, search array, resume) instead of one."
        ),
    )
    of_run.add_argument(
        "--script-dir",
        default="slurm_scripts",
        help=(
            "Directory to write generated scripts (used with --multinode, "
            "default: slurm_scripts)"
        ),
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "init":
        from convgeno.cli.init_cmd import run_init

        run_init(output_path=args.output)

    elif args.command == "orthofinder":
        if args.orthofinder_command == "generate":
            if args.multinode:
                from convgeno.cli.orthofinder_cmd import run_generate_multinode

                run_generate_multinode(
                    config_path=args.config, script_dir=args.script_dir
                )
            else:
                from convgeno.cli.orthofinder_cmd import run_generate

                run_generate(config_path=args.config, script_path=args.script)
        elif args.orthofinder_command == "run":
            if args.multinode:
                from convgeno.cli.orthofinder_cmd import run_submit_multinode

                run_submit_multinode(
                    config_path=args.config,
                    script_dir=args.script_dir,
                    skip_confirm=args.yes,
                )
            else:
                from convgeno.cli.orthofinder_cmd import run_submit

                run_submit(
                    config_path=args.config,
                    script_path=args.script,
                    skip_confirm=args.yes,
                )
        else:
            of_parser.print_help()
