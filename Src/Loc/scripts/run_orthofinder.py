"""
run_orthofinder.py — Run OrthoFinder on cleaned proteomes.

Inputs:  Data/interim/cleaned_proteomes/  (one FASTA per species)
Outputs: Data/processed/orthofinder/      (full OrthoFinder Results directory)
Logs:    Data/logs/orthofinder/

Copies cleaned proteomes into orthofinder_input/ (OrthoFinder requires
all FASTAs in a single directory), invokes OrthoFinder, then moves
results to processed/.

Usage:
    python Src/Loc/scripts/run_orthofinder.py \\
        --config Src/Loc/configs/example_config.yaml \\
        --tools  Src/Loc/configs/tool_paths.yaml
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Run OrthoFinder")
    parser.add_argument("--config", required=True)
    parser.add_argument("--tools", required=True)
    args = parser.parse_args()
    # TODO: implement — delegates to convgeno.external.orthofinder
    print("run_orthofinder.py: not yet implemented", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
