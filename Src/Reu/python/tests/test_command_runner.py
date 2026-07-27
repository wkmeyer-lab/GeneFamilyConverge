"""Tests for convgeno.utils.command_runner."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from convgeno.utils.command_runner import (
    CommandError,
    check_tool_available,
    run,
)

# Small helper scripts run via the current interpreter.  Written to files
# (not passed via ``-c``) to avoid shell-quoting pitfalls; only real newlines
# as separators, and ``print()`` so no backslash escapes are needed.
HELLO = "print('hello')\n"
FAIL = "import sys\nprint('boom', file=sys.stderr)\nsys.exit(3)\n"
ENVCHK = "import os\nprint(os.environ.get('CONVGENO_TEST', 'MISSING'))\n"


def _script(tmp_path: Path, name: str, body: str) -> Path:
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


class TestRun:
    def test_captures_stdout(self, tmp_path: Path):
        script = _script(tmp_path, "hello.py", HELLO)
        result = run([sys.executable, str(script)])
        assert result.returncode == 0
        assert result.stdout.strip() == "hello"

    def test_empty_command_raises(self):
        with pytest.raises(ValueError, match="non-empty command"):
            run([])

    def test_dry_run_does_not_execute(self, tmp_path: Path):
        marker = tmp_path / "marker.txt"
        body = f"from pathlib import Path\nPath(r'{marker}').write_text('x')\n"
        script = _script(tmp_path, "mk.py", body)

        result = run([sys.executable, str(script)], dry_run=True)

        assert result.returncode == 0
        assert result.stdout == ""
        assert not marker.exists()  # side effect never happened

    def test_nonzero_exit_raises_command_error(self, tmp_path: Path):
        script = _script(tmp_path, "fail.py", FAIL)
        with pytest.raises(CommandError, match="exit 3") as excinfo:
            run([sys.executable, str(script)])
        assert excinfo.value.returncode == 3
        assert "boom" in str(excinfo.value)  # stderr tail included

    def test_check_false_returns_instead_of_raising(self, tmp_path: Path):
        script = _script(tmp_path, "fail.py", FAIL)
        result = run([sys.executable, str(script)], check=False)
        assert result.returncode == 3

    def test_log_files_written_when_log_dir_given(self, tmp_path: Path):
        script = _script(tmp_path, "hello.py", HELLO)
        log_dir = tmp_path / "logs"
        run([sys.executable, str(script)], log_dir=log_dir)

        out_files = list(log_dir.glob("*.out"))
        err_files = list(log_dir.glob("*.err"))
        assert out_files and err_files
        assert "hello" in out_files[0].read_text(encoding="utf-8")

    def test_env_overrides_are_merged(self, tmp_path: Path):
        script = _script(tmp_path, "envchk.py", ENVCHK)
        result = run([sys.executable, str(script)], env={"CONVGENO_TEST": "xyz"})
        # Child saw our override *and* still ran (PATH etc. preserved).
        assert result.stdout.strip() == "xyz"


class TestCheckToolAvailable:
    def test_available_by_absolute_path(self):
        assert check_tool_available(sys.executable) is True

    def test_unavailable_returns_false(self):
        assert check_tool_available("nonexistent_tool_xyz_123") is False
