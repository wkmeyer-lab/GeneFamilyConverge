"""
validate_inputs.py — Verify all pipeline inputs before starting expensive runs.

Checks performed:
    - Every species in the phenotype table has a matching FASTA in proteome_dir
    - Every FASTA is parseable and non-empty
    - The species tree is valid Newick and its tip labels match the species set
    - The phenotype table has no missing values for required columns
    - Tool binaries are callable (validates tool_paths.yaml)

Usage:
    python Src/Loc/scripts/validate_inputs.py \\
        --config Src/Loc/configs/example_config.yaml \\
        --tools  Src/Loc/configs/tool_paths.yaml

Exit codes:
    0  All checks passed
    1  One or more checks failed (details printed to stderr)
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Validate pipeline inputs")
    parser.add_argument("--config", required=True, help="Path to pipeline config YAML")
    parser.add_argument("--tools", required=True, help="Path to tool_paths.yaml")
    args = parser.parse_args()
    # TODO: implement
    print("validate_inputs.py: not yet implemented", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
