"""
convgeno.external.orthofinder — Reusable wrapper for OrthoFinder CLI.

Functions to implement:
    build_command(proteome_dir, threads, extra_args, tool_path) -> list[str]
        Construct the OrthoFinder command as a list of args.
        Does NOT run it — that's command_runner's job.

    validate_output(output_dir) -> list[str]
        Check that expected OrthoFinder output files exist.
        Returns error messages for anything missing.

    find_results_dir(output_dir) -> Path
        OrthoFinder creates Results_<date>/ inside its output dir.
        This finds the most recent one.

    parse_gene_counts(results_dir) -> dict
        Read Orthogroups_GeneCount.tsv from a Results directory.
        Delegates actual parsing to convgeno.io.tables.
"""
