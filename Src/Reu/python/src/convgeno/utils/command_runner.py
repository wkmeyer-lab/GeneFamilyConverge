"""
convgeno.utils.command_runner — Single point of entry for all subprocess calls.

Every external tool invocation in the pipeline goes through this module.
No other file should import subprocess directly.

Functions to implement:
    run(cmd, log_dir=None, dry_run=False) -> subprocess.CompletedProcess
        Execute a command (list of strings).
        - Logs the full command before execution
        - Captures stdout and stderr
        - Writes stdout/stderr to log files if log_dir is provided
        - Raises a clear error on non-zero exit code
        - If dry_run=True, logs the command but does not execute it

    check_tool_available(tool_command) -> bool
        Check whether a tool binary is callable (e.g., `which orthofinder`).

Design notes:
    Centralizing subprocess calls here gives you:
    - One place to add dry-run support
    - One place to add timing/profiling
    - One place to log every command for reproducibility
    - One place to handle HPC-specific quirks (module load, etc.)
"""
