"""SLURM partition discovery and per-node CPU detection via sinfo."""

from __future__ import annotations

import os
import re
import stat
import subprocess
import tempfile
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


def _parse_slurm_memory_mb(value: str) -> int | None:
    raw = value.strip().rstrip("+")
    if raw.isdigit():
        parsed = int(raw)
        return parsed if parsed > 0 else None
    return None


def _parse_scontrol_memory_field(text: str, field: str) -> int | None:
    match = re.search(rf"\b{field}=([^\s]+)", text)
    if match is None:
        return None
    return _parse_slurm_memory_mb(match.group(1))


def detect_partition_memory(partition: str) -> dict[str, int | None]:
    """Detect memory limits and node memory for a SLURM partition.

    Returns ``None`` values when SLURM commands are unavailable or a field
    is absent. The function never raises, so ``convgeno init`` can fall back
    to an explicit user prompt.
    """
    detected: dict[str, int | None] = {
        "max_mem_per_cpu_mb": None,
        "def_mem_per_cpu_mb": None,
        "max_mem_per_node_mb": None,
        "def_mem_per_node_mb": None,
        "min_node_memory_mb": None,
    }

    try:
        partition_result = subprocess.run(
            ["scontrol", "show", "partition", partition],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if partition_result.returncode == 0:
            text = partition_result.stdout
            detected["max_mem_per_cpu_mb"] = _parse_scontrol_memory_field(
                text, "MaxMemPerCPU"
            )
            detected["def_mem_per_cpu_mb"] = _parse_scontrol_memory_field(
                text, "DefMemPerCPU"
            )
            detected["max_mem_per_node_mb"] = _parse_scontrol_memory_field(
                text, "MaxMemPerNode"
            )
            detected["def_mem_per_node_mb"] = _parse_scontrol_memory_field(
                text, "DefMemPerNode"
            )
    except Exception:
        pass

    try:
        node_result = subprocess.run(
            ["sinfo", "-N", "-h", "-p", partition, "-o", "%N %m"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if node_result.returncode == 0:
            node_memory: list[int] = []
            for line in node_result.stdout.splitlines():
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                memory_mb = _parse_slurm_memory_mb(parts[1])
                if memory_mb is not None:
                    node_memory.append(memory_mb)
            if node_memory:
                detected["min_node_memory_mb"] = min(node_memory)
    except Exception:
        pass

    return detected


def recommend_memory_mb(
    cpus_per_task: int,
    max_mem_per_cpu_mb: int | None,
    def_mem_per_cpu_mb: int | None,
    min_node_memory_mb: int | None,
) -> int:
    """Recommend an explicit total SLURM memory request in megabytes."""
    if max_mem_per_cpu_mb is not None and max_mem_per_cpu_mb > 0:
        return cpus_per_task * max_mem_per_cpu_mb

    if def_mem_per_cpu_mb is not None and def_mem_per_cpu_mb > 0:
        return cpus_per_task * def_mem_per_cpu_mb

    if min_node_memory_mb is not None and min_node_memory_mb > 0:
        return int(min_node_memory_mb * 0.90)

    raise ValueError(
        "Could not detect partition memory limits. Ask the user for an "
        "explicit SLURM memory request."
    )


def _empty_scratch_detection() -> dict[str, str | bool | None]:
    return {
        "scratch_base": None,
        "is_ephemeral": False,
        "method": "none",
        "write_granted": False,
    }


def _is_ephemeral_scratch_path(path: str) -> bool:
    return path.startswith(("/local", "/tmp", "/dev/shm"))


def _dir_is_writable(path: str) -> bool:
    """Definitive writability test: create then remove a probe file in ``path``.

    More reliable than ``os.access`` over NFS/root-squash filesystems, where
    ``os.access`` can report a permission the kernel will actually deny.
    """
    try:
        with tempfile.NamedTemporaryFile(dir=path, prefix=".convgeno_wtest_"):
            pass
        return True
    except OSError:
        return False


def _current_user_owns(path: str) -> bool:
    """Return ``True`` if ``path`` is owned by the effective user.

    Guarded for platforms without ``os.geteuid`` (e.g. Windows dev machines),
    where ownership cannot be determined and we conservatively return ``False``.
    """
    try:
        return os.stat(path).st_uid == os.geteuid()
    except (AttributeError, OSError):
        return False


def _grant_owner_write(path: str) -> bool:
    """Add the owner-write bit (equivalent to ``chmod u+w``) to ``path``.

    Returns ``True`` on success. Only the owner-write bit is added; group and
    other permissions are left untouched.
    """
    try:
        os.chmod(path, os.stat(path).st_mode | stat.S_IWUSR)
        return True
    except OSError:
        return False


def _usable_scratch_path(
    path: str, *, create_if_possible: bool
) -> tuple[str | None, bool]:
    """Return ``(usable_path, write_granted)`` for a scratch candidate.

    ``usable_path`` is ``path`` when it is (or can be made) a writable
    directory, else ``None``. ``write_granted`` is ``True`` only when an
    owner-write bit had to be added to make an owned-but-unwritable directory
    usable.
    """
    try:
        if not os.path.isdir(path):
            if not create_if_possible:
                return None, False
            parent = os.path.dirname(path) or "/"
            if not (os.path.isdir(parent) and os.access(parent, os.W_OK)):
                return None, False
            os.makedirs(path, exist_ok=True)
            if not os.path.isdir(path):
                return None, False

        # Directory exists; decide writability by a real write probe.
        if _dir_is_writable(path):
            return path, False

        # Exists but not writable — try to grant owner-write if we own it.
        if (
            _current_user_owns(path)
            and _grant_owner_write(path)
            and _dir_is_writable(path)
        ):
            return path, True

        return None, False
    except Exception:
        return None, False


def _scratch_path_is_on_home_filesystem(path: str) -> bool:
    try:
        return os.stat(path).st_dev == os.stat("/home").st_dev
    except Exception:
        return False


def detect_scratch_dir() -> dict[str, str | bool | None]:
    """Detect a writable scratch space root for cluster jobs.

    Environment-provided scratch is preferred, then common shared and
    node-local scratch locations are probed. Filesystem errors are ignored
    so init can continue on non-SLURM or restricted systems.
    """
    # ## NEW: Honor an explicit environment-provided scratch directory first.
    scratch_env = os.environ.get("SCRATCH")
    if scratch_env:
        scratch_base, write_granted = _usable_scratch_path(
            scratch_env, create_if_possible=False
        )
        if scratch_base is not None:
            return {
                "scratch_base": scratch_base,
                "is_ephemeral": _is_ephemeral_scratch_path(scratch_base),
                "method": "env_var",
                "write_granted": write_granted,
            }

    username = os.environ.get("USER", "unknown")
    candidates = [
        (f"/share/ceph/scratch/{username}", False),
        (f"/scratch/{username}", None),
        ("/local/scratch", True),
        ("/tmp/scratch", True),
    ]

    # ## NEW: Probe known cluster scratch paths without raising on failures.
    for path, fixed_ephemeral in candidates:
        scratch_base, write_granted = _usable_scratch_path(
            path, create_if_possible=True
        )
        if scratch_base is None:
            continue
        is_ephemeral = (
            not _scratch_path_is_on_home_filesystem(scratch_base)
            if fixed_ephemeral is None
            else fixed_ephemeral
        )
        return {
            "scratch_base": scratch_base,
            "is_ephemeral": is_ephemeral,
            "method": "path_probe",
            "write_granted": write_granted,
        }

    return _empty_scratch_detection()
