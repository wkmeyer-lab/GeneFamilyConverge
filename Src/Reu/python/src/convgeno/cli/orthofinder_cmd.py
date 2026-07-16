"""CLI handlers for ``convgeno orthofinder generate`` and ``convgeno orthofinder run``."""

from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from convgeno.slurm.config import PipelineConfig, normalize_optional_account


@dataclass
class SubmitResult:
    """Outcome of a single ``submit_sbatch`` call."""

    job_id: str
    account_stripped: bool = False
from convgeno.slurm.multinode_generator import (
    derive_search_array_sizing,
    generate_prepare_script,
    generate_resume_script,
    generate_search_array_script,
)
from convgeno.slurm.script_generator import generate_orthofinder_script, write_script
from convgeno.validation.orthofinder_inputs import validate_orthofinder_inputs


def validate_orthofinder_output_dir(output_dir: str) -> None:
    """Validate that ``output_dir`` is suitable for a fresh OrthoFinder run.

    OrthoFinder refuses to run when the non-default ``-o`` output
    directory already exists. We also reject paths containing literal
    quote characters, which usually indicate a malformed config value.
    The parent directory is created if missing.

    Raises
    ------
    ValueError
        If ``output_dir`` contains quote characters, or if the directory
        itself already exists.
    """
    if "'" in output_dir or '"' in output_dir:
        raise ValueError(
            f"Invalid output_dir {output_dir!r}: path must not contain "
            "quote characters."
        )

    path = Path(output_dir)
    if path.exists():
        raise ValueError(
            f"OrthoFinder output_dir already exists: {path}\n"
            "Choose a fresh output_dir. OrthoFinder refuses to use an "
            "existing non-default output directory."
        )

    path.parent.mkdir(parents=True, exist_ok=True)


def is_invalid_account_error(stderr: str) -> bool:
    """Return ``True`` if *stderr* indicates an invalid SLURM account."""
    text = stderr.lower()
    return (
        "invalid account" in text
        or "account/partition combination" in text
    )


_ACCOUNT_RE = re.compile(r"^\s*#SBATCH\s+--account=.*$", re.MULTILINE)


def strip_account_directive(script_text: str) -> str:
    """Return *script_text* with any ``#SBATCH --account=...`` line removed."""
    return _ACCOUNT_RE.sub("", script_text)


def _parse_job_id(stdout: str) -> str:
    """Extract the numeric job ID from ``sbatch`` stdout.

    Raises ``RuntimeError`` if the ID cannot be parsed.
    """
    job_id = stdout.strip().split()[-1] if stdout.strip() else ""
    if not job_id or not job_id.isdigit():
        raise RuntimeError(
            f"Could not parse job ID from sbatch output: {stdout.strip()}"
        )
    return job_id


def _run_sbatch(
    script_path: Path | str,
    dependency: str | None = None,
) -> subprocess.CompletedProcess:
    """Low-level sbatch invocation — no retry or error interpretation."""
    cmd = ["sbatch"]
    if dependency is not None:
        cmd.append(f"--dependency=afterok:{dependency}")
    cmd.append(str(script_path))
    return subprocess.run(cmd, capture_output=True, text=True, timeout=30)


def _submit_with_stripped_account(
    path: Path,
    dependency: str | None,
) -> SubmitResult:
    """Submit *path* after removing the ``#SBATCH --account`` line."""
    original_text = path.read_text(encoding="utf-8")
    stripped = strip_account_directive(original_text)
    if stripped != original_text:
        tmp = Path(tempfile.mktemp(
            suffix=".sh",
            dir=path.parent,
            prefix=f"{path.stem}_noaccount_",
        ))
        tmp.write_text(stripped, encoding="utf-8")
        tmp.chmod(tmp.stat().st_mode | 0o755)
        result = _run_sbatch(tmp, dependency)
        if result.returncode != 0:
            raise RuntimeError(
                f"sbatch submission failed (exit code {result.returncode}):\n"
                f"{result.stderr.strip()}"
            )
        return SubmitResult(job_id=_parse_job_id(result.stdout), account_stripped=True)
    # No account line to strip — submit as-is.
    result = _run_sbatch(path, dependency)
    if result.returncode != 0:
        raise RuntimeError(
            f"sbatch submission failed (exit code {result.returncode}):\n"
            f"{result.stderr.strip()}"
        )
    return SubmitResult(job_id=_parse_job_id(result.stdout), account_stripped=False)


def submit_sbatch(
    script_path: Path,
    dependency: str | None = None,
    *,
    strip_account: bool = False,
) -> SubmitResult:
    """Submit a SLURM batch script via sbatch with automatic account fallback.

    If the script contains a ``#SBATCH --account=…`` line and SLURM
    rejects the account, the submission is retried once with the account
    line stripped.  When *strip_account* is ``True`` the account line is
    removed up-front (used by multi-node callers after the first
    fallback succeeded).

    Returns a :class:`SubmitResult` with the job ID and whether the
    account directive was stripped during fallback.
    """
    path = Path(script_path)

    if strip_account:
        return _submit_with_stripped_account(path, dependency)

    result = _run_sbatch(path, dependency)

    if result.returncode == 0:
        return SubmitResult(job_id=_parse_job_id(result.stdout))

    script_text = path.read_text(encoding="utf-8")
    has_account_line = bool(_ACCOUNT_RE.search(script_text))

    # Case A: script has an account line and SLURM rejected it.
    if is_invalid_account_error(result.stderr) and has_account_line:
        account_match = _ACCOUNT_RE.search(script_text)
        account_val = ""
        if account_match:
            account_val = account_match.group().split("=", 1)[-1].strip()
        print(
            f"WARNING: SLURM rejected account '{account_val}' for this partition.\n"
            "Retrying submission without --account because some clusters "
            "do not require account names."
        )
        stripped = strip_account_directive(script_text)
        tmp = Path(tempfile.mktemp(
            suffix=".sh",
            dir=path.parent,
            prefix=f"{path.stem}_noaccount_",
        ))
        tmp.write_text(stripped, encoding="utf-8")
        tmp.chmod(tmp.stat().st_mode | 0o755)

        retry = _run_sbatch(tmp, dependency)
        if retry.returncode == 0:
            print(
                "Submitted successfully without a SLURM account.\n"
                "To avoid this warning in future runs, set "
                "slurm.account: null in pipeline_config.yaml."
            )
            return SubmitResult(
                job_id=_parse_job_id(retry.stdout),
                account_stripped=True,
            )

        raise RuntimeError(
            "Job submission failed with and without the configured account.\n"
            "This cluster may require a valid SLURM allocation account.\n"
            "Please provide a valid account name in pipeline_config.yaml "
            "under slurm.account.\n"
            f"Original account attempted: {account_val}\n"
            f"Original stderr:\n{result.stderr.strip()}\n"
            f"Retry stderr:\n{retry.stderr.strip()}"
        )

    # Case B: no account line but SLURM hints one is required.
    if not has_account_line and "account" in result.stderr.lower():
        raise RuntimeError(
            "This cluster appears to require a SLURM allocation account.\n"
            "Set a valid account in pipeline_config.yaml:\n\n"
            "slurm:\n"
            "  account: <your_valid_allocation>\n\n"
            f"sbatch stderr:\n{result.stderr.strip()}"
        )

    raise RuntimeError(
        f"sbatch submission failed (exit code {result.returncode}):\n"
        f"{result.stderr.strip()}"
    )


def run_generate(config_path: str, script_path: str) -> Path:
    """Load config, generate the OrthoFinder SLURM script, and write it to disk.

    Returns the script path.
    """
    config = PipelineConfig.load(config_path)

    if config.runtime is None:
        print(
            "ERROR: No runtime configuration found in config file.\n"
            "Run 'convgeno init' to detect and store conda paths."
        )
        sys.exit(1)

    if config.orthofinder is not None:
        validation = validate_orthofinder_inputs(config.orthofinder.input_dir)
        print(validation.summary())
        if not validation.is_valid():
            print("Fix the errors above before generating the SLURM script.")
            sys.exit(1)

        try:
            validate_orthofinder_output_dir(config.orthofinder.output_dir)
        except ValueError as exc:
            print(f"ERROR: {exc}")
            sys.exit(1)

    content = generate_orthofinder_script(config, config.runtime)
    path = write_script(content, script_path)
    print(f"SLURM script written to: {path}")
    return path


def run_submit(
    config_path: str,
    script_path: str,
    skip_confirm: bool = False,
) -> None:
    """Generate the OrthoFinder SLURM script, optionally show it for approval, and submit via sbatch."""
    path = run_generate(config_path, script_path)

    if not skip_confirm:
        print(f"\n--- Generated script: {path} ---")
        print(path.read_text(encoding="utf-8"))
        print("--- End of script ---\n")
        response = input("Submit this job? [Y/n]: ").strip().lower()
        if response.startswith("n"):
            print(f"Aborted. Script saved to: {path}")
            return

    try:
        sr = submit_sbatch(path)
    except FileNotFoundError:
        print("ERROR: sbatch command not found. Are you on a SLURM login node?")
        sys.exit(1)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("ERROR: sbatch submission timed out.")
        sys.exit(1)

    print("Job submitted successfully.")
    print(f"  Job ID:  {sr.job_id}")
    print(f"  Script:  {path}")
    print(f"  Monitor: squeue -j {sr.job_id}")
    print(f"  Cancel:  scancel {sr.job_id}")


def run_generate_multinode(config_path: str, script_dir: str) -> dict:
    """Load config, validate inputs, and generate the multi-node scripts.

    Returns a dict mapping script role to the on-disk script path. The v2
    search array is not yet implemented, so only ``"prepare"`` and ``"resume"``
    are generated; the search phase is reported as pending.
    """
    config = PipelineConfig.load(config_path)

    if config.runtime is None:
        print(
            "ERROR: No runtime configuration found in config file.\n"
            "Run 'convgeno init' to detect and store conda paths."
        )
        sys.exit(1)

    if config.orthofinder is not None:
        validation = validate_orthofinder_inputs(config.orthofinder.input_dir)
        print(validation.summary())
        if not validation.is_valid():
            print("Fix the errors above before generating the SLURM scripts.")
            sys.exit(1)

        try:
            validate_orthofinder_output_dir(config.orthofinder.output_dir)
        except ValueError as exc:
            print(f"ERROR: {exc}")
            sys.exit(1)

    # Derive the search-array sizing ONCE and share it with both generators so
    # the manifest bucket count (prepare) and the --array width (search) agree.
    sizing = derive_search_array_sizing(config)
    prepare_content = generate_prepare_script(config, config.runtime, sizing=sizing)
    search_content = generate_search_array_script(config, config.runtime, sizing=sizing)
    resume_content = generate_resume_script(config, config.runtime)

    script_dir_path = Path(script_dir)
    prepare_path = write_script(prepare_content, script_dir_path / "orthofinder_prepare.sh")
    search_path = write_script(search_content, script_dir_path / "orthofinder_search.sh")
    resume_path = write_script(resume_content, script_dir_path / "orthofinder_resume.sh")

    throttle = min(sizing.concurrency, sizing.tasks)
    print(f"Multi-node SLURM scripts written to {script_dir}/")
    print(f"  Prepare: {prepare_path}")
    print(
        f"  Search:  {search_path}  "
        f"(--array=0-{sizing.tasks - 1}%{throttle}, {sizing.cpus} cpus/task, "
        f"--time {sizing.time_limit}, --mem {sizing.mem})"
    )
    print(f"  Resume:  {resume_path}")

    return {
        "prepare": prepare_path,
        "search": search_path,
        "resume": resume_path,
    }


def run_submit_multinode(
    config_path: str,
    script_dir: str,
    skip_confirm: bool = False,
) -> None:
    """Generate and submit the three-job multi-node OrthoFinder dependency chain."""
    paths = run_generate_multinode(config_path, script_dir)

    if not skip_confirm:
        for role in ("prepare", "search", "resume"):
            print(f"\n--- {role.title()} script: {paths[role]} ---")
            print(paths[role].read_text(encoding="utf-8"))
            print(f"--- End of {role} script ---\n")
        response = input("Submit this 3-job dependency chain? [Y/n]: ").strip().lower()
        if response.startswith("n"):
            print(f"Aborted. Scripts saved in: {script_dir}")
            return

    try:
        # If the first submission triggers an invalid-account fallback and
        # succeeds without an account, submit the remaining scripts without
        # the account line as well so the dependency chain is consistent.
        sr_a = submit_sbatch(paths["prepare"])
        strip_rest = sr_a.account_stripped
        print(f"Job A (prepare) submitted: {sr_a.job_id}")

        sr_b = submit_sbatch(
            paths["search"],
            dependency=sr_a.job_id,
            strip_account=strip_rest,
        )
        if sr_b.account_stripped:
            strip_rest = True
        print(f"Job B (search array) submitted: {sr_b.job_id}")

        sr_c = submit_sbatch(
            paths["resume"],
            dependency=sr_b.job_id,
            strip_account=strip_rest,
        )
        print(f"Job C (resume) submitted: {sr_c.job_id}")
    except FileNotFoundError:
        print("ERROR: sbatch command not found. Are you on a SLURM login node?")
        sys.exit(1)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("ERROR: sbatch submission timed out.")
        sys.exit(1)

    print("")
    print("Multi-node OrthoFinder job chain submitted:")
    print(f"  Job A (prepare):      {sr_a.job_id}")
    print(f"  Job B (search array): {sr_b.job_id} (depends on {sr_a.job_id})")
    print(f"  Job C (resume):       {sr_c.job_id} (depends on {sr_b.job_id})")
    print("")
    print("Monitor all: squeue -u $USER")
    print(f"Cancel all:  scancel {sr_a.job_id} {sr_b.job_id} {sr_c.job_id}")
