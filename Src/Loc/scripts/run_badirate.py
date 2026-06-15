"""
run_badirate.py — Run BadiRate to compare null vs phenotype-driven models.

For each gene family, fits null and alternative evolution models, then
compares via AIC to identify families whose evolution is linked to the
phenotype.

Inputs:  Data/interim/badirate_input/
Outputs: Data/processed/badirate/
Logs:    Data/logs/badirate/

Usage:
    python Src/Loc/scripts/run_badirate.py \\
        --config Src/Loc/configs/example_config.yaml \\
        --tools  Src/Loc/configs/tool_paths.yaml
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Run BadiRate")
    parser.add_argument("--config", required=True)
    parser.add_argument("--tools", required=True)
    args = parser.parse_args()
    # TODO: implement — delegates to convgeno.external.badirate
    print("run_badirate.py: not yet implemented", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
