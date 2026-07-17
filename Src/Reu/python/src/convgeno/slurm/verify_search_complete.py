"""Completeness gate for the multi-node search phase.

The resume SLURM script invokes this BEFORE ``orthofinder -b``. It enumerates
the ``n²`` expected ``Blast{i}_{j}.txt.gz`` results from ``SpeciesIDs.txt`` and
checks that every one exists (and is non-empty) in the WorkingDirectory. If any
are missing it prints them, maps them back to the search-array task ids that
would produce them (by scanning the per-task manifests), prints a ready-to-paste
``sbatch --array=<ids>`` resubmit line, and exits non-zero WITHOUT touching the
results -- so ``-b`` never runs against an incomplete search set. Run as::

    python -m convgeno.slurm.verify_search_complete \\
        --work-dir <WorkingDirectory> \\
        --manifest-dir <RUN>_search_manifests \\
        [--search-script orthofinder_search.sh]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from convgeno.slurm.multinode_generator import (
    _parse_species_pair,
    read_species_ids,
)

_MANIFEST_NAME_RE = re.compile(r"search_task_(\d+)\.txt$")
_MAX_LISTED = 50  # cap how many missing pairs are printed in full


def find_missing_blast_files(
    work_dir: Path,
    ids: list[int],
) -> list[tuple[int, int]]:
    """Return the ``(i, j)`` pairs whose ``Blast{i}_{j}.txt.gz`` is missing.

    "Missing" means the file is absent or empty (zero bytes). Deeper gzip
    integrity is already enforced upstream: the search array skips a result
    only when it passes ``gzip -t``, and the ``afterok`` dependency lets resume
    start only when every search task succeeded.
    """
    work_dir = Path(work_dir)
    missing: list[tuple[int, int]] = []
    for i in ids:
        for j in ids:
            f = work_dir / f"Blast{i}_{j}.txt.gz"
            if not (f.is_file() and f.stat().st_size > 0):
                missing.append((i, j))
    return missing


def map_pairs_to_task_ids(
    manifest_dir: Path,
    pairs: list[tuple[int, int]],
) -> list[int]:
    """Array task ids whose manifests contain any of ``pairs`` (ascending).

    The species pair of each manifest command is parsed to integers (via the
    ``Blast{i}_{j}`` token), so matching is exact -- no ``Blast0_1`` vs
    ``Blast0_10`` substring collision.
    """
    wanted = set(pairs)
    manifest_dir = Path(manifest_dir)
    tasks: set[int] = set()
    if not manifest_dir.is_dir():
        return []
    for manifest in manifest_dir.glob("search_task_*.txt"):
        name_match = _MANIFEST_NAME_RE.search(manifest.name)
        if name_match is None:
            continue
        task_id = int(name_match.group(1))
        for line in manifest.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                pair = _parse_species_pair(line)
            except ValueError:
                continue
            if pair in wanted:
                tasks.add(task_id)
                break
    return sorted(tasks)


def verify(
    work_dir: Path,
    manifest_dir: Path,
    search_script: str,
) -> tuple[bool, str]:
    """Check completeness. Returns ``(ok, report)``.

    ``ok`` is ``True`` only when all ``n²`` results are present.
    """
    ids = read_species_ids(work_dir)
    expected = len(ids) * len(ids)
    missing = find_missing_blast_files(work_dir, ids)
    present = expected - len(missing)

    lines = [
        f"Completeness gate: {present}/{expected} Blast results present "
        f"({len(ids)} species)."
    ]
    if not missing:
        lines.append("All search results present; proceeding to resume.")
        return True, "\n".join(lines)

    lines.append(f"ERROR: {len(missing)} of {expected} search results missing.")
    for i, j in missing[:_MAX_LISTED]:
        lines.append(f"  missing Blast{i}_{j}.txt.gz")
    if len(missing) > _MAX_LISTED:
        lines.append(f"  ... and {len(missing) - _MAX_LISTED} more")

    task_ids = map_pairs_to_task_ids(manifest_dir, missing)
    array_spec = ",".join(str(t) for t in task_ids) if task_ids else "<task ids>"
    lines.append("")
    lines.append(
        "Resubmit the affected search-array tasks, then rerun this resume job:"
    )
    lines.append(f"  sbatch --array={array_spec} {search_script}")
    lines.append("Results were NOT touched; orthofinder -b did not run.")
    return False, "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point. Returns 0 when complete, 1 otherwise."""
    parser = argparse.ArgumentParser(
        prog="python -m convgeno.slurm.verify_search_complete",
        description="Verify all n^2 search results exist before resume (-b).",
    )
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--manifest-dir", required=True)
    parser.add_argument("--search-script", default="orthofinder_search.sh")
    args = parser.parse_args(argv)

    try:
        ok, report = verify(
            Path(args.work_dir),
            Path(args.manifest_dir),
            args.search_script,
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    stream = sys.stdout if ok else sys.stderr
    print(report, file=stream)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
