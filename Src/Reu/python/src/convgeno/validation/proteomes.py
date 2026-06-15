"""
convgeno.validation.proteomes — Validate proteome FASTA inputs.

Functions to implement:
    validate_proteome(path, source_format) -> list[str]
        Check a single proteome FASTA for problems.
        Returns a list of error/warning messages (empty = valid).
        Checks: file is parseable, non-empty, headers match expected format,
        no duplicate sequence IDs.

    validate_proteome_set(directory, species_list, source_format) -> dict[str, list[str]]
        Validate all proteomes for a species set.
        Returns {species_name: [errors]} for any species with problems.
        Also checks that every species in the list has a matching file.
"""
