"""Entry point for the convgeno command-line interface."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _resolve_orthofinder_mode(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> str:
    """Resolve the effective OrthoFinder execution mode.

    Default is ``multinode``. Accepts ``--mode {multinode,single-node}`` plus
    the back-compat aliases ``--multinode`` / ``--single-node``, and errors on
    a conflicting combination.
    """
    wants_single = args.single_node or args.mode == "single-node"
    wants_multi = args.multinode or args.mode == "multinode"
    if wants_single and wants_multi:
        parser.error(
            "Conflicting execution mode: choose only one of "
            "--mode / --multinode / --single-node."
        )
    return "single-node" if wants_single else "multinode"


def _resolve_config_path(explicit: str | None, mode: str) -> str:
    """Config path to use: explicit ``--config`` else the mode-specific default.

    ``multinode`` -> ``pipeline_config_multinode.yaml``; ``single-node`` ->
    ``pipeline_config_singlenode.yaml`` (both written by ``convgeno init``).
    Falls back to the legacy ``pipeline_config.yaml`` when the mode-specific
    file is absent.
    """
    if explicit is not None:
        return explicit
    mode_default = (
        "pipeline_config_singlenode.yaml"
        if mode == "single-node"
        else "pipeline_config_multinode.yaml"
    )
    return mode_default if Path(mode_default).exists() else "pipeline_config.yaml"


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
        default=None,
        help=(
            "Path to pipeline config file (default: the mode-specific "
            "pipeline_config_<mode>.yaml, else pipeline_config.yaml)"
        ),
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
        "--mode",
        choices=["multinode", "single-node"],
        default=None,
        help=(
            "Execution mode (default: multinode — prepare + search array + "
            "resume). 'single-node' generates the single benchmark/fallback "
            "script."
        ),
    )
    of_gen.add_argument(
        "--multinode",
        action="store_true",
        help="Alias for --mode multinode (this is the default).",
    )
    of_gen.add_argument(
        "--single-node",
        action="store_true",
        help="Alias for --mode single-node (benchmark/fallback).",
    )
    of_gen.add_argument(
        "--script-dir",
        default="slurm_scripts",
        help=(
            "Directory to write multi-node scripts (default: slurm_scripts)"
        ),
    )

    of_run = of_subparsers.add_parser(
        "run",
        help="Generate the OrthoFinder SLURM script and submit it.",
    )
    of_run.add_argument(
        "--config",
        default=None,
        help=(
            "Path to pipeline config file (default: the mode-specific "
            "pipeline_config_<mode>.yaml, else pipeline_config.yaml)"
        ),
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
        "--mode",
        choices=["multinode", "single-node"],
        default=None,
        help=(
            "Execution mode (default: multinode — submit prepare + search "
            "array + resume). 'single-node' submits the single "
            "benchmark/fallback job."
        ),
    )
    of_run.add_argument(
        "--multinode",
        action="store_true",
        help="Alias for --mode multinode (this is the default).",
    )
    of_run.add_argument(
        "--single-node",
        action="store_true",
        help="Alias for --mode single-node (benchmark/fallback).",
    )
    of_run.add_argument(
        "--script-dir",
        default="slurm_scripts",
        help=(
            "Directory to write multi-node scripts (default: slurm_scripts)"
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
            mode = _resolve_orthofinder_mode(parser, args)
            config_path = _resolve_config_path(args.config, mode)
            if mode == "single-node":
                from convgeno.cli.orthofinder_cmd import run_generate

                run_generate(config_path=config_path, script_path=args.script)
            else:
                from convgeno.cli.orthofinder_cmd import run_generate_multinode

                run_generate_multinode(
                    config_path=config_path, script_dir=args.script_dir
                )
        elif args.orthofinder_command == "run":
            mode = _resolve_orthofinder_mode(parser, args)
            config_path = _resolve_config_path(args.config, mode)
            if mode == "single-node":
                from convgeno.cli.orthofinder_cmd import run_submit

                run_submit(
                    config_path=config_path,
                    script_path=args.script,
                    skip_confirm=args.yes,
                )
            else:
                from convgeno.cli.orthofinder_cmd import run_submit_multinode

                run_submit_multinode(
                    config_path=config_path,
                    script_dir=args.script_dir,
                    skip_confirm=args.yes,
                )
        else:
            of_parser.print_help()
