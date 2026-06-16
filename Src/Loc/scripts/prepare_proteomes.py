"""
prepare_proteomes.py — Filter proteomes to longest isoform per gene.

For each species FASTA in proteome_dir:
    1. Parse headers to extract gene ID and transcript ID
    2. Group sequences by gene ID
    3. Keep only the longest protein sequence per gene
    4. Write filtered FASTA to interim/cleaned_proteomes/

The filtered proteomes become the input to OrthoFinder.

Usage:
    python Src/Loc/scripts/prepare_proteomes.py \\
        --config Src/Loc/configs/example_config.yaml
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Filter proteomes to longest isoform")
    parser.add_argument("--config", required=True, help="Path to pipeline config YAML")
    args = parser.parse_args()
    # TODO: implement — delegates to convgeno.io.fasta and convgeno.validation.proteomes
    print("prepare_proteomes.py: not yet implemented", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
