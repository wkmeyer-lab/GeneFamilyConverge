"""Build the balanced per-task search manifest at prepare time.

The prepare SLURM script invokes this AFTER ``orthofinder -op`` has produced the
WorkingDirectory (``Species*.fa``) and the blastp-only ``_search_commands.txt``.
It LPT-partitions the ``n^2`` commands into ``T`` cost-balanced buckets (``T``
fixed at generate time so the array width matches) and writes one manifest file
per array task. Run as::

    python -m convgeno.slurm.build_search_manifest \\
        --search-commands <RUN>_search_commands.txt \\
        --work-dir <WorkingDirectory> \\
        --manifest-dir <RUN>_search_manifests \\
        --tasks <T>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from convgeno.slurm.multinode_generator import (
    build_search_commands,
    lpt_partition,
    read_species_fasta_sizes,
    write_task_manifests,
)


def build_manifest(
    search_commands: Path,
    work_dir: Path,
    manifest_dir: Path,
    tasks: int,
) -> dict:
    """Partition the emitted search commands into ``tasks`` manifests.

    Returns a summary dict (command/task counts + the LPT sizing byproducts).
    The same ``command_lines`` list is passed to ``build_search_commands`` and
    ``write_task_manifests`` so bucket line indices align exactly.
    """
    command_lines = search_commands.read_text(encoding="utf-8").splitlines()
    if not any(line.strip() for line in command_lines):
        raise ValueError(f"No search commands found in {search_commands}")

    species_sizes = read_species_fasta_sizes(work_dir)
    commands = build_search_commands(command_lines, species_sizes)
    result = lpt_partition(commands, tasks)
    paths = write_task_manifests(result.buckets, command_lines, manifest_dir)

    return {
        "commands": len(commands),
        "tasks": len(paths),
        "max_bucket_cost": result.max_bucket_cost,
        "mem_determinant_bytes": result.mem_determinant_bytes,
        "manifest_dir": str(manifest_dir),
    }


def main(argv: list[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code (0 on success)."""
    parser = argparse.ArgumentParser(
        prog="python -m convgeno.slurm.build_search_manifest",
        description="Build the balanced per-task search manifest (LPT).",
    )
    parser.add_argument("--search-commands", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--manifest-dir", required=True)
    parser.add_argument("--tasks", type=int, required=True)
    args = parser.parse_args(argv)

    try:
        summary = build_manifest(
            Path(args.search_commands),
            Path(args.work_dir),
            Path(args.manifest_dir),
            args.tasks,
        )
    except (ValueError, FileNotFoundError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"Wrote {summary['tasks']} task manifests for "
        f"{summary['commands']} search commands to {summary['manifest_dir']}"
    )
    print(
        f"  biggest bucket cost = {summary['max_bucket_cost']} bytes^2; "
        f"largest DB = {summary['mem_determinant_bytes']} bytes"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
