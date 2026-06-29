"""SLURM partition discovery and per-node CPU detection via sinfo."""

from __future__ import annotations

import re
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


def _empty_node_cpu_detection() -> dict[str, int]:
    return {
        "min_cpus_per_node": 0,
        "max_cpus_per_node": 0,
        "recommended_cpus": 0,
        "threads_per_core": 1,
        "physical_cores": 0,
        "recommended_physical": 16,
        "node_count": 0,
    }


def detect_node_cpus(partition: str) -> dict[str, int]:
    """Detect CPU counts and physical-core defaults for a SLURM partition.

    The recommendation uses the smallest node in the partition so jobs can
    land on any node. It reserves four physical cores as memory-bandwidth
    headroom, with a floor of eight cores, and falls back to the historical
    static default when SLURM commands are unavailable.
    """
    try:
        sinfo_result = subprocess.run(
            ["sinfo", "-p", partition, "-N", "--noheader", "-o", "%N %c"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception:
        return _empty_node_cpu_detection()

    if sinfo_result.returncode != 0:
        return _empty_node_cpu_detection()

    node_names: list[str] = []
    cpu_counts: list[int] = []
    for line in sinfo_result.stdout.splitlines():
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        try:
            cpu_count = int(parts[1])
        except ValueError:
            continue
        node_names.append(parts[0])
        cpu_counts.append(cpu_count)

    if not cpu_counts:
        return _empty_node_cpu_detection()

    threads_per_core = 1
    try:
        scontrol_result = subprocess.run(
            ["scontrol", "show", "node", node_names[0]],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if scontrol_result.returncode == 0:
            match = re.search(r"\bThreadsPerCore=(\d+)\b", scontrol_result.stdout)
            if match is not None:
                detected_threads = int(match.group(1))
                if detected_threads > 0:
                    threads_per_core = detected_threads
    except Exception:
        threads_per_core = 1

    min_cpus = min(cpu_counts)
    max_cpus = max(cpu_counts)
    physical_cores = min_cpus // threads_per_core
    return {
        "min_cpus_per_node": min_cpus,
        "max_cpus_per_node": max_cpus,
        "recommended_cpus": max(min_cpus - 4, 8),
        "threads_per_core": threads_per_core,
        "physical_cores": physical_cores,
        "recommended_physical": max(physical_cores - 4, 8),
        "node_count": len(cpu_counts),
    }
