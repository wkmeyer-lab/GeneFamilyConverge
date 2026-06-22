"""Configuration dataclasses for external bioinformatics tools."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field


@dataclass(frozen=True)
class OrthoFinderConfig:
    """Configuration for OrthoFinder gene family clustering."""

    input_dir: str
    output_dir: str
    search_threads: int = 16
    analysis_threads: int = 8
    sequence_search: str = "diamond"
    msa_program: str = "mafft"
    tree_program: str = "fasttree"
    extra_args: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Serialize all fields to a plain dict."""
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> OrthoFinderConfig:
        """Construct from a plain dict (e.g. parsed from YAML).

        Unknown keys are silently ignored so config files remain
        forward-compatible across package versions.

        Raises
        ------
        ValueError
            If required keys ``input_dir`` or ``output_dir`` are missing.
        """
        if "input_dir" not in data:
            raise ValueError("OrthoFinderConfig requires 'input_dir' to be set.")
        if "output_dir" not in data:
            raise ValueError("OrthoFinderConfig requires 'output_dir' to be set.")
        valid_fields = {f.name for f in dataclasses.fields(cls)}
        extra_args = data.pop("extra_args", [])
        filtered = {k: v for k, v in data.items() if k in valid_fields}
        return cls(extra_args=extra_args, **filtered)

    def to_command_args(self) -> list[str]:
        """Build the command-line argument list for the ``orthofinder`` command."""
        args = [
            "-f", self.input_dir,
            "-o", self.output_dir,
            "-t", str(self.search_threads),
            "-a", str(self.analysis_threads),
            "-S", self.sequence_search,
            "-M", self.msa_program,
            "-T", self.tree_program,
        ]
        args.extend(self.extra_args)
        return args
