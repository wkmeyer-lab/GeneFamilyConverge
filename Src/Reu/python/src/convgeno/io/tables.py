"""
convgeno.io.tables — Read and write TSV/CSV tables consistently.

Functions to implement:
    read_tsv(path) -> list[dict]
        Read a TSV into a list of row dicts.

    write_tsv(rows, path, columns)
        Write a list of row dicts to a TSV with specified column order.

    read_gene_counts(path) -> dict
        Read OrthoFinder's Orthogroups_GeneCount.tsv into a structured dict.

Design notes:
    - Thin wrappers around csv.DictReader/DictWriter.
    - Keeps pandas out of the core library to minimize dependencies.
      Loc scripts can use pandas if they want; Reu stays lightweight.
"""
