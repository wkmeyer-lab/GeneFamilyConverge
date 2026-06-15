"""
convgeno.utils.config — Load and validate YAML config files.

Functions to implement:
    load_config(path) -> dict
        Read a YAML config file and return its contents as a dict.

    load_tool_paths(path) -> dict
        Read tool_paths.yaml and return {tool_name: command_string}.

    resolve_paths(config, base_dir) -> dict
        Resolve all relative paths in a config dict against a base directory.
        Ensures paths work regardless of where the script is invoked from.

    validate_config(config) -> list[str]
        Check that required keys exist in the config.
        Returns list of error messages.
"""
