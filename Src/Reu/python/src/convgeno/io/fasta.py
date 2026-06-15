"""
convgeno.io.fasta — Read, write, and inspect FASTA files.

Functions to implement:
    parse_fasta(path) -> dict[str, str]
        Read a FASTA file into {header: sequence} dict.

    write_fasta(records, path)
        Write {header: sequence} dict to a FASTA file.

    extract_gene_id(header, source_format) -> str
        Extract the gene ID from a FASTA header.
        source_format selects the parsing logic (e.g., "ensembl").

    extract_transcript_id(header, source_format) -> str
        Extract the transcript/protein ID from a FASTA header.

    filter_longest_isoform(records, source_format) -> dict[str, str]
        Group by gene ID, keep longest sequence per gene.

Design notes:
    - Uses BioPython SeqIO under the hood.
    - source_format dispatches to a parser function; only "ensembl" is
      implemented initially. Add others when needed.
    - All functions accept paths (str or Path), not open file handles,
      to keep the interface simple.
"""
