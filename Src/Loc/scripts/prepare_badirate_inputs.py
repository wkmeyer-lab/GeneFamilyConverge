"""
prepare_badirate_inputs.py — Convert OrthoFinder output + phenotype data to BadiRate format.

BadiRate compares two models of gene family evolution:
    - Null:  uniform turnover rate across the whole tree
    - Alt:   separate turnover rates per phenotype state

This script prepares:
    - Gene family size data per species (from OrthoFinder)
    - A phylogenetic tree with branch lengths
    - A phenotype-labeled tree for the alternative model

Inputs:  Data/processed/orthofinder/
         config phenotype table
Outputs: Data/interim/badirate_input/

Usage:
    python Src/Loc/scripts/prepare_badirate_inputs.py \\
        --config Src/Loc/configs/example_config.yaml
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Prepare BadiRate inputs")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    # TODO: implement
    print("prepare_badirate_inputs.py: not yet implemented", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
