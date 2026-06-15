"""
convgeno.external.cafe — Reusable wrapper for CAFE-5 CLI.

Functions to implement:
    build_command(count_file, tree_file, n_gamma_cats, extra_args, tool_path) -> list[str]
        Construct the CAFE-5 command.

    validate_input(count_file, tree_file) -> list[str]
        Check that CAFE-5 input files exist and are correctly formatted.
        CAFE-5 expects: tab-delimited with Desc and Family ID columns,
        all counts as non-negative integers, and an ultrametric Newick tree.

    validate_output(output_dir) -> list[str]
        Check expected CAFE-5 output files exist.

    format_gene_counts_for_cafe(orthofinder_counts_path, output_path)
        Reformat OrthoFinder gene count table to CAFE-5's expected layout.
        Adds Desc column, adjusts column names.
"""
