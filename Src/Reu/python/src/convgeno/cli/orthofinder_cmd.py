"""CLI handlers for ``convgeno orthofinder generate`` and ``convgeno orthofinder run``."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from convgeno.slurm.config import PipelineConfig
from convgeno.slurm.multinode_generator import (
    generate_prepare_script,
    generate_resume_script,
    generate_search_array_script,
)
from convgeno.slurm.script_generator import generate_orthofinder_script, write_script
from convgeno.validation.orthofinder_inputs import validate_orthofinder_inputs


def submit_sbatch(script_path: Path, dependency: str | None = None) -> str:
    """Submit a SLURM batch script via sbatch.

    Optionally adds an ``afterok`` dependency on a previously-submitted
    job ID. Returns the new job ID as a string.

    Raises
    ------
    RuntimeError
        If submission fails or the job ID cannot be parsed from the output.
    """
    cmd = ["sbatch"]
    if dependency is not None:
        cmd.append(f"--dependency=afterok:{dependency}")
    cmd.append(str(script_path))
    result = subprocess.run(
        cmd,
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


def run_generate_multinode(config_path: str, script_dir: str) -> dict:
    """Load config, validate inputs, and generate the three multi-node scripts.

    Returns a dict mapping script role (``"prepare"``, ``"search"``,
    ``"resume"``) to the on-disk script path.
    """
    config = PipelineConfig.load(config_path)

    if config.orthofinder is not None:
        validation = validate_orthofinder_inputs(config.orthofinder.input_dir)
        print(validation.summary())
        if not validation.is_valid():
            print("Fix the errors above before generating the SLURM scripts.")
            sys.exit(1)

    prepare_content = generate_prepare_script(config)
    search_content = generate_search_array_script(config)
    resume_content = generate_resume_script(config)

    script_dir_path = Path(script_dir)
    prepare_path = write_script(prepare_content, script_dir_path / "orthofinder_prepare.sh")
    search_path = write_script(search_content, script_dir_path / "orthofinder_search_array.sh")
    resume_path = write_script(resume_content, script_dir_path / "orthofinder_resume.sh")

    print(f"Multi-node SLURM scripts written to {script_dir}/")
    print(f"  Prepare: {prepare_path}")
    print(f"  Search:  {search_path}")
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
        job_a_id = submit_sbatch(paths["prepare"])
        print(f"Job A (prepare) submitted: {job_a_id}")
        job_b_id = submit_sbatch(paths["search"], dependency=job_a_id)
        print(f"Job B (search array) submitted: {job_b_id}")
        job_c_id = submit_sbatch(paths["resume"], dependency=job_b_id)
        print(f"Job C (resume) submitted: {job_c_id}")
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
    print(f"  Job A (prepare):      {job_a_id}")
    print(f"  Job B (search array): {job_b_id} (depends on {job_a_id})")
    print(f"  Job C (resume):       {job_c_id} (depends on {job_b_id})")
    print("")
    print("Monitor all: squeue -u $USER")
    print(f"Cancel all:  scancel {job_a_id} {job_b_id} {job_c_id}")
