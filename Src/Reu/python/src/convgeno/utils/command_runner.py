"""
convgeno.utils.command_runner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Single point of entry for all subprocess calls in the pipeline.

Every external tool invocation (OrthoFinder, r8s, CAFE-5, …) goes through
:func:`run`.  No other module should import :mod:`subprocess` directly, so that
there is exactly one place to add dry-run support, timing, per-command logging
for reproducibility, and HPC-specific quirks (module load, etc.).

Public API
----------
- :func:`run` — execute a command, capture output, optionally write log files,
  and raise a clear :class:`CommandError` on failure.
- :func:`check_tool_available` — is a tool binary on ``PATH`` (or a valid
  executable path)?
- :class:`CommandError` — raised for a non-zero exit status.
"""

from __future__ import annotations

import logging
import os
import re
import shlex
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

logger = logging.getLogger(__name__)

__all__ = ["CommandError", "run", "check_tool_available"]

#: Max stderr lines appended to a :class:`CommandError` message.
_STDERR_TAIL_LINES = 20


class CommandError(RuntimeError):
    """Raised when an external command exits with a non-zero status.

    Attributes
    ----------
    cmd : list[str]
        The command that was executed.
    returncode : int
        The process exit status.
    stdout, stderr : str or None
        Captured output (if any).
    """

    def __init__(
        self,
        cmd: Sequence[str],
        returncode: int,
        stdout: str | None = None,
        stderr: str | None = None,
    ) -> None:
        self.cmd = [str(part) for part in cmd]
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr

        message = f"Command failed (exit {returncode}): {_display(self.cmd)}"
        if stderr and stderr.strip():
            tail = stderr.strip().splitlines()[-_STDERR_TAIL_LINES:]
            message += "\n--- stderr (tail) ---\n" + "\n".join(tail)
        super().__init__(message)


def _display(cmd: Sequence[str]) -> str:
    """Render a command list as a copy-pasteable shell string for logs."""
    return shlex.join(str(part) for part in cmd)


def _write_logs(
    log_dir: Path, cmd: Sequence[str], completed: subprocess.CompletedProcess
) -> None:
    """Write captured stdout/stderr to ``<log_dir>/<tool>.out`` / ``.err``."""
    log_dir.mkdir(parents=True, exist_ok=True)
    tool = Path(str(cmd[0])).name or "command"
    stem = re.sub(r"[^A-Za-z0-9._-]", "_", tool)
    (log_dir / f"{stem}.out").write_text(completed.stdout or "", encoding="utf-8")
    (log_dir / f"{stem}.err").write_text(completed.stderr or "", encoding="utf-8")


def run(
    cmd: Sequence[str],
    log_dir: Path | str | None = None,
    dry_run: bool = False,
    *,
    cwd: Path | str | None = None,
    env: Mapping[str, str] | None = None,
    timeout: float | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Execute *cmd* and return the completed process.

    The command is logged before execution.  stdout and stderr are always
    captured (as text) on the returned :class:`subprocess.CompletedProcess`.

    Parameters
    ----------
    cmd : sequence of str
        The command and its arguments.  Never run through a shell
        (``shell=False``), so no quoting/escaping is required or performed.
    log_dir : Path or str, optional
        If given, stdout/stderr are also written to ``<log_dir>/<tool>.out``
        and ``.err`` (the directory is created if needed).
    dry_run : bool
        If True, log the command and return an empty successful
        ``CompletedProcess`` without executing anything.
    cwd : Path or str, optional
        Working directory for the child process.
    env : mapping, optional
        Environment overrides, merged on top of the current environment
        (so ``PATH`` etc. are preserved).
    timeout : float, optional
        Seconds before the child is killed and ``TimeoutExpired`` is raised.
    check : bool
        If True (default), raise :class:`CommandError` on a non-zero exit.

    Returns
    -------
    subprocess.CompletedProcess
        With ``.returncode``, ``.stdout``, ``.stderr`` (text).

    Raises
    ------
    ValueError
        If *cmd* is empty.
    CommandError
        If the command exits non-zero and ``check`` is True.
    subprocess.TimeoutExpired
        If *timeout* is exceeded.
    """
    cmd = [str(part) for part in cmd]
    if not cmd:
        raise ValueError("run() requires a non-empty command.")

    display = _display(cmd)
    logger.info("Running: %s", display)

    if dry_run:
        logger.info("Dry run: command not executed.")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    run_env = {**os.environ, **env} if env is not None else None

    completed = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(cwd) if cwd is not None else None,
        env=run_env,
        timeout=timeout,
        check=False,
    )

    if log_dir is not None:
        _write_logs(Path(log_dir), cmd, completed)

    if completed.returncode != 0:
        logger.error("Command exited with status %d: %s", completed.returncode, display)
        if check:
            raise CommandError(
                cmd, completed.returncode, completed.stdout, completed.stderr
            )
    else:
        logger.debug("Command succeeded: %s", display)

    return completed


def check_tool_available(tool_command: str) -> bool:
    """Return True if *tool_command* is callable.

    Uses :func:`shutil.which`, so it accepts either a bare command name resolved
    against ``PATH`` (e.g. ``"orthofinder"``) or an explicit path to an
    executable (e.g. ``"/opt/r8s/r8s"``).

    Parameters
    ----------
    tool_command : str
        Command name or path to check.

    Returns
    -------
    bool
        Whether the tool was found and is executable.
    """
    resolved = shutil.which(str(tool_command))
    if resolved is not None:
        logger.debug("Tool available: %s -> %s", tool_command, resolved)
        return True
    logger.debug("Tool NOT found on PATH: %s", tool_command)
    return False
