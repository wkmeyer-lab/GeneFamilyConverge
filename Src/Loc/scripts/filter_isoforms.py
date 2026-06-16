#!/usr/bin/env python3
"""
Step 1 · Filter proteome FASTA files to the longest isoform per gene.

Two modes of operation
----------------------
**Directory mode** (batch — one file per species)::

    python step01_filter_isoforms.py dir data/raw/proteomes/ data/interim/primary_transcripts/

**File mode** (single species, one or more input files)::

    python step01_filter_isoforms.py files \\
        -i data/raw/proteomes/Homo_sapiens.pep.all.fa.gz \\
        -o data/interim/primary_transcripts/Homo_sapiens.fa

File mode also handles per-chromosome downloads for a single species::

    python step01_filter_isoforms.py files \\
        -i proteome_chr1.fa.gz proteome_chr2.fa.gz \\
        -o data/interim/primary_transcripts/species.fa

Snakemake rules should call **file mode** (one species per rule
invocation); directory mode is for interactive / manual use.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from convgeno.io.fasta import filter_longest_isoforms, process_directory
from convgeno.utils.logging import setup_logging

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------
# Argument parsing
# -------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="step01_filter_isoforms",
        description="Filter proteome FASTA files to the longest isoform per gene.",
    )

    # --- shared options ---
    parser.add_argument(
        "--format",
        choices=["auto", "ensembl", "ncbi"],
        default="auto",
        help="Header format for gene-ID extraction (default: auto-detect).",
    )
    parser.add_argument(
        "--on-duplicate",
        choices=["error", "warn", "skip"],
        default="error",
        help="How to handle duplicate accession lines (default: error).",
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

    sub = parser.add_subparsers(dest="mode", required=True)

    # --- directory mode ---
    dir_p = sub.add_parser(
        "dir",
        help="Batch: process every FASTA file in a directory.",
    )
    dir_p.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing one FASTA file per species.",
    )
    dir_p.add_argument(
        "output_dir",
        type=Path,
        help="Directory for filtered output files.",
    )
    dir_p.add_argument(
        "--stats-json",
        type=Path,
        default=None,
        help="Write per-species summary statistics to this JSON file.",
    )

    # --- file mode ---
    file_p = sub.add_parser(
        "files",
        help="Single species: process one or more FASTA files.",
    )
    file_p.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=Path,
        required=True,
        dest="input_files",
        help="One or more input FASTA files (all for the same species).",
    )
    file_p.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output FASTA file path.",
    )

    return parser


# -------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run the isoform filter.  Returns exit code."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    setup_logging(level=args.log_level, log_file=args.log_file)

    try:
        if args.mode == "dir":
            all_stats = process_directory(
                input_dir=args.input_dir,
                output_dir=args.output_dir,
                header_format=args.format,
                on_duplicate=args.on_duplicate,
            )
            # Print summary table.
            print(f"\n{'File':<50} {'Total':>7} {'Genes':>7} {'Unid.':>7} {'Written':>7}")
            print("-" * 80)
            for fname, st in sorted(all_stats.items()):
                print(
                    f"{fname:<50} "
                    f"{st['total_records']:>7} "
                    f"{st['unique_genes']:>7} "
                    f"{st['unidentified']:>7} "
                    f"{st['output_records']:>7}"
                )
            # Optionally dump full stats to JSON.
            if args.stats_json:
                args.stats_json.parent.mkdir(parents=True, exist_ok=True)
                with open(args.stats_json, "w") as jf:
                    json.dump(all_stats, jf, indent=2)
                logger.info("Wrote stats to %s", args.stats_json)

        elif args.mode == "files":
            stats = filter_longest_isoforms(
                input_paths=args.input_files,
                output_path=args.output,
                header_format=args.format,
                on_duplicate=args.on_duplicate,
            )
            print(
                f"Done: {stats['total_records']} records → "
                f"{stats['unique_genes']} genes + "
                f"{stats['unidentified']} unidentified = "
                f"{stats['output_records']} written to {args.output}"
            )

    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        logger.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())