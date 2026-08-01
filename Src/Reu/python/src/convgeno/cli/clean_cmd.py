"""``convgeno clean``: the pipeline's first step — longest-isoform filtering.

Filters a directory of raw proteomes (one FASTA per species) to the longest
isoform per gene, writing the cleaned set into the directory OrthoFinder runs
on. Inputs come either from explicit ``raw_dir``/``out_dir`` arguments or, when
those are omitted, from the ``proteome_input`` block of ``pipeline_config.yaml``
(written by ``convgeno init``). This is what the Snakemake ``clean_proteomes``
rule invokes, and it doubles as a standalone command.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from convgeno.io.fasta import process_directory
from convgeno.utils.logging import setup_logging

_CONFIG_CANDIDATES = (
    "pipeline_config_multinode.yaml",
    "pipeline_config_singlenode.yaml",
    "pipeline_config.yaml",
)


def _default_config_path() -> str | None:
    """First existing default config, or ``None`` if none is present.

    ``proteome_input`` is identical across the mode-specific configs, so any of
    them serves for cleaning.
    """
    for candidate in _CONFIG_CANDIDATES:
        if Path(candidate).exists():
            return candidate
    return None


def run_clean(
    config_path: str | None = None,
    raw_dir: str | None = None,
    out_dir: str | None = None,
    header_format: str | None = None,
    on_duplicate: str | None = None,
    stats_json: str | None = None,
    log_level: str = "INFO",
) -> int:
    """Run longest-isoform filtering. Returns a process exit code."""
    setup_logging(level=log_level)

    # When raw/out are not both given, fill the gaps from the config's
    # proteome_input block.
    if raw_dir is None or out_dir is None:
        path = config_path or _default_config_path()
        if path is None:
            print(
                "No raw_dir/out_dir given and no pipeline config found. "
                "Run 'convgeno init' first, or pass RAW_DIR and OUT_DIR.",
                file=sys.stderr,
            )
            return 1
        from convgeno.slurm.config import PipelineConfig

        try:
            config = PipelineConfig.load(path)
        except (FileNotFoundError, ValueError) as exc:
            print(f"Could not read config {path}: {exc}", file=sys.stderr)
            return 1
        pi = config.proteome_input
        if pi is None or not pi.has_raw():
            print(
                f"Config {path} has no 'proteome_input.raw_dir' to clean. "
                "Re-run 'convgeno init' to set one, or pass RAW_DIR and OUT_DIR.",
                file=sys.stderr,
            )
            return 1
        raw_dir = raw_dir or pi.raw_dir
        out_dir = out_dir or pi.cleaned_dir
        header_format = header_format or pi.header_format
        on_duplicate = on_duplicate or pi.on_duplicate

    header_format = header_format or "auto"
    on_duplicate = on_duplicate or "error"

    try:
        all_stats = process_directory(
            input_dir=Path(raw_dir),
            output_dir=Path(out_dir),
            header_format=header_format,
            on_duplicate=on_duplicate,
        )
    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        print(f"Cleaning failed: {exc}", file=sys.stderr)
        return 1

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
    print(f"\nCleaned {len(all_stats)} proteome(s) -> {out_dir}")

    if stats_json:
        stats_path = Path(stats_json)
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        with open(stats_path, "w", encoding="utf-8") as jf:
            json.dump(all_stats, jf, indent=2)
        print(f"Wrote per-species stats to {stats_path}")

    return 0
