"""Tests for convgeno.validation.proteomes"""

# Priority test cases:
#   - Valid proteome -> empty error list
#   - Empty file -> error reported
#   - Duplicate sequence IDs -> error reported
#   - Species in list but no matching FASTA -> error reported
#   - FASTA exists but species not in list -> warning reported
