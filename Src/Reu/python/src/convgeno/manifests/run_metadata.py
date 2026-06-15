"""
convgeno.manifests.run_metadata — Record pipeline run metadata.

Functions to implement:
    record_run_metadata(config_path, tool_paths, output_path)
        Write a YAML file capturing:
        - Timestamp
        - Config file used (and its checksum)
        - Tool versions (by calling each tool's --version flag)
        - Python version
        - R version
        - Hostname
        - Git commit hash (if in a repo)

    record_tool_version(tool_command) -> str
        Call a tool with --version / -v / --help and extract version string.
"""
