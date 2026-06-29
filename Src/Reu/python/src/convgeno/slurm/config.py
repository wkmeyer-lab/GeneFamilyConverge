"""SLURM and pipeline configuration dataclasses with YAML serialization."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.runtime import (
    CondaRuntimeConfig,
    runtime_config_from_dict,
    runtime_config_to_dict,
)


def normalize_optional_account(account: str | None) -> str | None:
    """Normalize an optional SLURM account value.

    Returns ``None`` for blank, ``"null"``, ``"None"``, or actual
    ``None`` inputs so that no ``#SBATCH --account`` line is emitted.
    A real non-empty string is returned as-is.
    """
    if account is None:
        return None
    value = str(account).strip()
    if value == "" or value.lower() in {"null", "none"}:
        return None
    return value


@dataclass(frozen=True)
class SlurmConfig:
    """SLURM resource configuration for a single job submission."""

    partition: str
    time_limit: str = "48:00:00"
    nodes: int = 1
    ntasks: int = 1
    cpus_per_task: int = 16
    mem_per_cpu: str = "4G"
    account: Optional[str] = None
    mail_user: Optional[str] = None
    mail_type: str = "END,FAIL"
    output_pattern: str = "logs/%x_%j.out"
    error_pattern: str = "logs/%x_%j.err"
    open_file_limit: Optional[int] = None
    extra_sbatch_args: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Serialize to a dict, omitting fields whose value is ``None``."""
        return {k: v for k, v in dataclasses.asdict(self).items() if v is not None}

    @classmethod
    def from_dict(cls, data: dict) -> SlurmConfig:
        """Construct from a plain dict (e.g. parsed from YAML).

        Unknown keys are silently ignored so config files remain
        forward-compatible across package versions.

        Raises
        ------
        ValueError
            If the required ``partition`` key is missing.
        """
        if "partition" not in data:
            raise ValueError(
                "SlurmConfig requires 'partition' to be set. "
                "Run 'convgeno init' to configure."
            )
        valid_fields = {f.name for f in dataclasses.fields(cls)}
        extra_sbatch_args = data.pop("extra_sbatch_args", [])
        filtered = {k: v for k, v in data.items() if k in valid_fields}
        if "account" in filtered:
            filtered["account"] = normalize_optional_account(filtered["account"])
        return cls(extra_sbatch_args=extra_sbatch_args, **filtered)

    def to_sbatch_lines(self) -> list[str]:
        """Generate ``#SBATCH`` directive lines for an sbatch script header.

        The account value is normalized through
        :func:`normalize_optional_account` so that ``"null"``, ``"None"``,
        and blank strings are treated identically to ``None`` — no
        ``--account`` line is emitted.
        """
        flag_map = {
            "partition": "--partition",
            "time_limit": "--time",
            "nodes": "--nodes",
            "ntasks": "--ntasks",
            "cpus_per_task": "--cpus-per-task",
            "mem_per_cpu": "--mem-per-cpu",
            "mail_user": "--mail-user",
            "mail_type": "--mail-type",
            "output_pattern": "--output",
            "error_pattern": "--error",
        }
        lines: list[str] = []
        for attr, flag in flag_map.items():
            value = getattr(self, attr)
            if value is not None:
                lines.append(f"#SBATCH {flag}={value}")

        account = normalize_optional_account(self.account)
        if account is not None:
            lines.append(f"#SBATCH --account={account}")

        for arg in self.extra_sbatch_args:
            lines.append(f"#SBATCH {arg}")
        return lines


@dataclass(frozen=True)
class PipelineConfig:
    """Top-level pipeline configuration combining project paths and SLURM settings."""

    project_dir: str
    conda_env: str = "convgeno"
    slurm: SlurmConfig = field(default_factory=lambda: SlurmConfig(partition=""))
    orthofinder: Optional[OrthoFinderConfig] = None
    runtime: Optional[CondaRuntimeConfig] = None

    def save(self, path: Path | str) -> None:
        """Write the configuration to a human-readable YAML file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        data: dict = {
            "project_dir": self.project_dir,
            "conda_env": self.conda_env,
            "slurm": self.slurm.to_dict(),
        }
        if self.orthofinder is not None:
            data["orthofinder"] = self.orthofinder.to_dict()
        if self.runtime is not None:
            data.update(runtime_config_to_dict(self.runtime))
        with open(path, "w", encoding="utf-8") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    @classmethod
    def load(cls, path: Path | str) -> PipelineConfig:
        """Load configuration from a YAML file.

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If required keys are missing from the YAML data.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(
                f"Config file not found: {path}. Run 'convgeno init' to create one."
            )
        with open(path, encoding="utf-8") as f:
            data = yaml.safe_load(f)
        if "project_dir" not in data:
            raise ValueError("Config file is missing required key 'project_dir'.")
        slurm_dict = data.get("slurm", {})
        slurm_config = SlurmConfig.from_dict(slurm_dict)
        of_dict = data.get("orthofinder")
        of_config = OrthoFinderConfig.from_dict(of_dict) if of_dict else None
        runtime_dict = data.get("runtime")
        runtime_config = (
            runtime_config_from_dict(runtime_dict) if runtime_dict else None
        )
        return cls(
            project_dir=data["project_dir"],
            conda_env=data.get("conda_env", "convgeno"),
            slurm=slurm_config,
            orthofinder=of_config,
            runtime=runtime_config,
        )
