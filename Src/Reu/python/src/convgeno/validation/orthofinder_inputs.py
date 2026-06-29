"""Pre-submission validation of OrthoFinder input directories."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

FASTA_EXTENSIONS = {".fa", ".fasta", ".faa"}
MIN_SPECIES_COUNT = 4
MAX_FILENAME_LENGTH = 200
FILENAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


@dataclass
class ValidationResult:
    """Collects errors and warnings from input validation.

    Errors are fatal — the pipeline should not proceed.
    Warnings are informational — the pipeline can proceed but the user
    should be aware.
    """

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    fasta_files: list[Path] = field(default_factory=list)

    def is_valid(self) -> bool:
        """Return ``True`` if no fatal errors were recorded."""
        return len(self.errors) == 0

    def summary(self) -> str:
        """Return a human-readable summary of the validation outcome."""
        lines: list[str] = []
        if self.errors:
            lines.append("VALIDATION FAILED")
            for msg in self.errors:
                lines.append(f"  ERROR: {msg}")
        else:
            lines.append("Validation passed.")
        for msg in self.warnings:
            lines.append(f"  WARNING: {msg}")
        lines.append(f"  FASTA files found: {len(self.fasta_files)}")
        return "\n".join(lines)


def validate_orthofinder_inputs(input_dir: str | Path) -> ValidationResult:
    """Validate the OrthoFinder input directory before script generation.

    Checks that the directory exists, contains enough FASTA files, that
    filenames are sane, and that files are non-empty and look like FASTA.
    Returns a ``ValidationResult`` with all errors and warnings collected.
    """
    input_dir = Path(input_dir).resolve()
    result = ValidationResult()

    if not input_dir.exists():
        result.errors.append(f"Input directory does not exist: {input_dir}")
        return result

    if not input_dir.is_dir():
        result.errors.append(f"Input path is not a directory: {input_dir}")
        return result

    fasta_files = sorted(
        f for f in input_dir.iterdir()
        if f.is_file() and f.suffix.lower() in FASTA_EXTENSIONS
    )
    result.fasta_files = fasta_files

    if not fasta_files:
        result.errors.append(
            f"No FASTA files found in {input_dir}. "
            f"Expected files with extensions: {', '.join(sorted(FASTA_EXTENSIONS))}"
        )
        return result

    if len(fasta_files) < MIN_SPECIES_COUNT:
        result.errors.append(
            f"Found {len(fasta_files)} FASTA file(s) but OrthoFinder requires "
            f"at least {MIN_SPECIES_COUNT} species."
        )

    stem_counts: Counter[str] = Counter()

    for fasta in fasta_files:
        if len(fasta.name) > MAX_FILENAME_LENGTH:
            result.warnings.append(
                f"Filename exceeds {MAX_FILENAME_LENGTH} characters: {fasta.name}"
            )

        if not FILENAME_PATTERN.match(fasta.stem):
            result.warnings.append(
                f"Filename contains unusual characters that may cause problems: "
                f"{fasta.name}. Recommended: alphanumeric, dots, hyphens, "
                f"underscores only."
            )

        try:
            size = fasta.stat().st_size
        except OSError:
            result.warnings.append(
                f"Could not read file to verify FASTA format: {fasta.name}"
            )
            stem_counts[fasta.stem] += 1
            continue

        if size == 0:
            result.errors.append(f"FASTA file is empty (0 bytes): {fasta.name}")
            stem_counts[fasta.stem] += 1
            continue

        try:
            with open(fasta, encoding="utf-8") as fh:
                first_line = fh.readline()
            if not first_line.startswith(">"):
                result.errors.append(
                    f"File does not start with a FASTA header (>): {fasta.name}"
                )
        except (UnicodeDecodeError, OSError):
            result.warnings.append(
                f"Could not read file to verify FASTA format: {fasta.name}"
            )

        stem_counts[fasta.stem] += 1

    for stem, count in stem_counts.items():
        if count > 1:
            result.errors.append(
                f"Duplicate species name detected: {stem}. "
                f"Each species must have exactly one FASTA file."
            )

    return result
