"""
run_cafe.py — Run CAFE-5 to model gene family evolution.

Estimates gene family size and turnover rate on each branch of the
phylogeny. Identifies gene families with significantly rapid evolution.

Inputs:  Data/interim/cafe_input/cafe_input.tsv
         Data/interim/cafe_input/species_tree.nwk
Outputs: Data/processed/cafe/
Logs:    Data/logs/cafe/

Usage:
    python Src/Loc/scripts/run_cafe.py \\
        --config Src/Loc/configs/example_config.yaml \\
        --tools  Src/Loc/configs/tool_paths.yaml
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Run CAFE-5")
    parser.add_argument("--config", required=True)
    parser.add_argument("--tools", required=True)
    args = parser.parse_args()
    # TODO: implement — delegates to convgeno.external.cafe
    print("run_cafe.py: not yet implemented", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
