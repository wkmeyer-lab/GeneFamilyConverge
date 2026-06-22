"""CLI handlers for ``convgeno orthofinder generate`` and ``convgeno orthofinder run``."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from convgeno.slurm.config import PipelineConfig
from convgeno.slurm.script_generator import generate_orthofinder_script, write_script
from convgeno.validation.orthofinder_inputs import validate_orthofinder_inputs


def submit_sbatch(script_path: Path) -> str:
    """Submit a SLURM batch script via sbatch.

    Returns the job ID as a string.

    Raises
    ------
    RuntimeError
        If submission fails or the job ID cannot be parsed from the output.
    """
    result = subprocess.run(
        ["sbatch", str(script_path)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"sbatch submission failed (exit code {result.returncode}):\n"
            f"{result.stderr.strip()}"
        )
    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else ""
    if not job_id or not job_id.isdigit():
        raise RuntimeError(
            f"Could not parse job ID from sbatch output: {result.stdout.strip()}"
        )
    return job_id


def run_generate(config_path: str, script_path: str) -> Path:
    """Load config, generate the OrthoFinder SLURM script, and write it to disk.

    Returns the script path.
    """
    config = PipelineConfig.load(config_path)

    if config.orthofinder is not None:
        validation = validate_orthofinder_inputs(config.orthofinder.input_dir)
        print(validation.summary())
        if not validation.is_valid():
            print("Fix the errors above before generating the SLURM script.")
            sys.exit(1)

    content = generate_orthofinder_script(config)
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
        job_id = submit_sbatch(path)
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
    print(f"  Job ID:  {job_id}")
    print(f"  Script:  {path}")
    print(f"  Monitor: squeue -j {job_id}")
    print(f"  Cancel:  scancel {job_id}")
