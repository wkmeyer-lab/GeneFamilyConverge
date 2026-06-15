"""Tests for convgeno.utils.command_runner"""

# Priority test cases:
#   - Run a simple command (e.g., echo) -> captures stdout
#   - Non-zero exit code -> raises clear exception
#   - dry_run=True -> logs command but does not execute
#   - check_tool_available("echo") -> True
#   - check_tool_available("nonexistent_tool_xyz") -> False
#   - Log files written when log_dir is provided
