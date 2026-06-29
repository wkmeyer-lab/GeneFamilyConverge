"""Configuration dataclasses for external bioinformatics tools."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class OrthoFinderConfig:
    """Configuration for OrthoFinder gene family clustering."""

    input_dir: str
    output_dir: str
    search_threads: int = 16
    analysis_threads: Optional[int] = 8
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
        """Build the command-line argument list for the ``orthofinder`` command.

        Notes
        -----
        OrthoFinder's ``-M`` flag takes a *method* name (``msa`` or
        ``dendroblast``), not a program name. The MSA program is selected
        with ``-A`` and the tree program with ``-T``. When ``msa_program``
        is truthy we emit ``-M msa -A <msa_program>`` (plus ``-T`` if a
        tree program is set). When ``msa_program`` is empty we fall back
        to OrthoFinder's default DendroBLAST gene-tree method and omit
        ``-M``, ``-A``, and ``-T`` entirely.
        """
        effective_analysis = self.analysis_threads if self.analysis_threads is not None else 1
        args = [
            "-f", self.input_dir,
            "-o", self.output_dir,
            "-t", str(self.search_threads),
            "-a", str(effective_analysis),
            "-S", self.sequence_search,
        ]
        if self.msa_program:
            args.extend(["-M", "msa", "-A", self.msa_program])
            if self.tree_program:
                args.extend(["-T", self.tree_program])
        args.extend(self.extra_args)
        return args
