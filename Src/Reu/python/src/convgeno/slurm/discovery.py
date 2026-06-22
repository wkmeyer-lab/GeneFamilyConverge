"""SLURM partition discovery via sinfo."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from typing import List


@dataclass
class PartitionInfo:
    """Information about a single SLURM partition discovered via sinfo."""

    name: str
    is_default: bool = False
    time_limit: str = "unknown"
    max_cpus_per_node: int = 0
    max_mem_mb_per_node: int = 0

    def __str__(self) -> str:
        label = f"{self.name} (default)" if self.is_default else self.name
        return (
            f"{label} — time limit: {self.time_limit}, "
            f"cpus/node: {self.max_cpus_per_node}, "
            f"mem/node: {self.max_mem_mb_per_node} MB"
        )


def _parse_sinfo_line(line: str) -> PartitionInfo | None:
    """Parse one line of ``sinfo -o "%P %l %c %m" --noheader`` output.

    Returns ``None`` if the line cannot be parsed (too few fields, etc.).
    """
    parts = line.strip().split()
    if len(parts) < 4:
        return None

    raw_name = parts[0]
    is_default = raw_name.endswith("*")
    name = raw_name.rstrip("*") if is_default else raw_name

    time_limit = parts[1]

    try:
        max_cpus = int(parts[2])
    except ValueError:
        max_cpus = 0

    try:
        max_mem = int(parts[3].rstrip("+"))
    except ValueError:
        max_mem = 0

    return PartitionInfo(
        name=name,
        is_default=is_default,
        time_limit=time_limit,
        max_cpus_per_node=max_cpus,
        max_mem_mb_per_node=max_mem,
    )


def discover_partitions() -> List[PartitionInfo]:
    """Discover available SLURM partitions by calling sinfo.

    Returns an empty list if sinfo is not available or fails.
    Never raises exceptions.
    """
    try:
        result = subprocess.run(
            ["sinfo", "-o", "%P %l %c %m", "--noheader"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode != 0:
            return []

        best: dict[str, PartitionInfo] = {}
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            info = _parse_sinfo_line(line)
            if info is None:
                continue
            existing = best.get(info.name)
            if existing is None or info.max_cpus_per_node > existing.max_cpus_per_node:
                best[info.name] = info

        partitions = list(best.values())
        partitions.sort(key=lambda p: (not p.is_default, p.name))
        return partitions
    except Exception:
        return []
