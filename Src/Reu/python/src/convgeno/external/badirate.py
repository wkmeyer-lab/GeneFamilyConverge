"""
convgeno.external.badirate — Reusable wrapper for BadiRate CLI.

Functions to implement:
    build_command(input_dir, tree_file, extra_args, tool_path) -> list[str]
        Construct the BadiRate command.

    validate_input(input_dir, tree_file) -> list[str]
        Check BadiRate inputs exist and are correctly formatted.

    validate_output(output_dir) -> list[str]
        Check expected BadiRate output files exist.

    format_inputs_for_badirate(gene_counts, tree_path, phenotype_map, output_dir)
        Prepare BadiRate-specific input files:
        - Gene family size data per species
        - Phenotype-labeled tree for alternative model

Notes:
    BadiRate is a Perl script (BadiRate.pl), not a compiled binary.
    Installation is manual — not available via conda/bioconda.
    See docs in Src/Loc/configs/tool_paths.yaml for path configuration.
"""
