"""
prepare_cafe_inputs.py — Convert OrthoFinder output to CAFE-5 input format.

CAFE-5 requires:
    1. A tab-delimited gene count file:
           Desc  Family_ID  Species1  Species2  ...
           (null)  OG0000001  3  5  ...
    2. A rooted ultrametric species tree in Newick format

This script:
    - Reads Orthogroups_GeneCount.tsv from OrthoFinder output
    - Reformats column names and adds the Desc column CAFE-5 expects
    - Validates / copies the species tree
    - Writes both to Data/interim/cafe_input/

Usage:
    python Src/Loc/scripts/prepare_cafe_inputs.py \\
        --config Src/Loc/configs/example_config.yaml
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Prepare CAFE-5 inputs")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    # TODO: implement
    print("prepare_cafe_inputs.py: not yet implemented", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
