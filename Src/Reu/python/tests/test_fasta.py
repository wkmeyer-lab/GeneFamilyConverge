"""Tests for convgeno.io.fasta"""

# Priority test cases:
#   - Parse a simple multi-sequence FASTA -> correct dict
#   - Extract gene ID from Ensembl header
#   - Handle empty file gracefully
#   - Handle malformed headers (missing gene ID field)
#   - Round-trip: parse then write produces identical file
#   - filter_longest_isoform: gene with 3 isoforms -> keeps longest
#   - filter_longest_isoform: gene with 1 isoform -> kept as-is
#   - filter_longest_isoform: tie in length -> deterministic pick
