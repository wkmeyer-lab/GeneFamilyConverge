"""``convgeno run``: drive the whole pipeline as one Snakemake DAG.

The DAG (``workflow/Snakefile``) is: clean proteomes -> OrthoFinder (+ ultrametric
and phenotype trees) -> CAFE-5 inputs -> run CAFE-5. This command renders a
Snakemake SLURM profile from ``pipeline_config.yaml`` and then either

* **submits a small, long-walltime orchestrator job** (default) that runs
  Snakemake, which fans the DAG out across the cluster; or
* runs Snakemake **in the foreground** (``--local``) on the current node — for an
  interactive allocation or debugging; or
* prints the planned DAG without running anything (``-n/--dry-run``).

Login nodes kill long-running processes (the same limit that OOM-kills a big
conda solve), so the default submits an orchestrator job rather than running
Snakemake on the login node.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import yaml

from convgeno.slurm.config import PipelineConfig, normalize_optional_account
from convgeno.slurm.runtime import render_conda_bootstrap


def _mem_to_mb(mem: str | None) -> int | None:
    """Parse a SLURM memory string ('350400M', '120G', '8000') to MB."""
    if not mem:
        return None
    text = str(mem).strip()
    units = {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 * 1024}
    suffix = text[-1].upper()
    if suffix in units:
        return max(1, int(round(float(text[:-1]) * units[suffix])))
    return max(1, int(round(float(text))))  # bare number -> assume MB


def _time_to_minutes(time_limit: str | None) -> int | None:
    """Parse a SLURM time string ('72:00:00', '3-00:00:00', '90:00') to minutes."""
    if not time_limit:
        return None
    text = str(time_limit).strip()
    days = 0
    if "-" in text:
        day_str, text = text.split("-", 1)
        days = int(day_str)
    parts = [int(p) for p in text.split(":")]
    if len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:  # SLURM 'MM:SS'
        hours, minutes, seconds = 0, parts[0], parts[1]
    elif len(parts) == 1:  # bare minutes
        hours, minutes, seconds = 0, parts[0], 0
    else:
        raise ValueError(f"Unrecognized SLURM time: {time_limit!r}")
    total = days * 1440 + hours * 60 + minutes + (1 if seconds else 0)
    return max(1, total)


def render_slurm_profile(config: PipelineConfig, jobs: int) -> dict:
    """Build the Snakemake SLURM-executor profile from a pipeline config.

    Only the SLURM-submitted rules (CAFE-5) draw on these defaults; the light
    steps (clean, cafe_inputs, orthofinder-wrapper) are local rules.
    """
    slurm = config.slurm
    default_resources: dict = {}
    if slurm.partition:
        default_resources["slurm_partition"] = slurm.partition
    account = normalize_optional_account(slurm.account)
    if account is not None:
        default_resources["slurm_account"] = account
    mem_mb = _mem_to_mb(slurm.mem)
    if mem_mb is not None:
        default_resources["mem_mb"] = mem_mb
    runtime = _time_to_minutes(slurm.time_limit)
    if runtime is not None:
        default_resources["runtime"] = runtime

    profile: dict = {"executor": "slurm", "jobs": jobs}
    if default_resources:
        profile["default-resources"] = default_resources
    return profile


def _snakemake_common(snakefile: Path, config_abs: Path, mode: str) -> list[str]:
    return [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--config",
        f"pipeline_config={config_abs}",
        f"mode={mode}",
    ]


def _orchestrator_script(
    config: PipelineConfig,
    snakemake_cmd: list[str],
    project_dir: Path,
) -> str:
    """Render the small sbatch script that runs Snakemake on a compute node."""
    slurm = config.slurm
    account = normalize_optional_account(slurm.account)
    # The script runs on the Linux cluster, so always emit POSIX paths even when
    # `convgeno run` is invoked from a Windows checkout.
    pd = project_dir.as_posix()
    header = [
        "#!/bin/bash",
        "#SBATCH --job-name=convgeno_run",
    ]
    if slurm.partition:
        header.append(f"#SBATCH --partition={slurm.partition}")
    if account is not None:
        header.append(f"#SBATCH --account={account}")
    header += [
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=1",
        "#SBATCH --mem=8G",
        f"#SBATCH --time={slurm.time_limit}",
        f"#SBATCH --output={pd}/logs/convgeno_run_%j.out",
        f"#SBATCH --error={pd}/logs/convgeno_run_%j.err",
        "",
        "set -euo pipefail",
        "",
    ]
    body = []
    if config.runtime is not None:
        body.append(render_conda_bootstrap(config.runtime))
        body.append("")
    body.append("# r8s and CAFE-5 are external tools installed under ~/.local/bin")
    body.append('export PATH="$HOME/.local/bin:$PATH"')
    body.append("")
    body.append(f'cd "{pd}"')
    body.append("")
    body.append(" ".join(_shquote(part) for part in snakemake_cmd))
    return "\n".join(header + body) + "\n"


def _shquote(arg: str) -> str:
    """Minimal shell-quote: wrap in double quotes if it contains spaces."""
    return f'"{arg}"' if " " in arg else arg


def run(
    config_path: str,
    mode: str,
    *,
    local: bool = False,
    dry_run: bool = False,
    jobs: int = 8,
    skip_confirm: bool = False,
) -> None:
    """Render the profile and run/submit the Snakemake DAG."""
    config = PipelineConfig.load(config_path)
    project_dir = Path(config.project_dir)
    config_abs = Path(config_path).resolve()

    snakefile = project_dir / "workflow" / "Snakefile"
    if not snakefile.exists():
        print(f"ERROR: workflow file not found: {snakefile}")
        sys.exit(1)

    # Dry run: show the DAG with the local executor (no SLURM plugin needed),
    # regardless of --local, and never submit anything.
    if dry_run:
        cmd = _snakemake_common(snakefile, config_abs, mode) + [
            "-n",
            "--cores",
            "1",
            "all",
        ]
        print("Dry run (no jobs will be submitted):")
        print("  " + " ".join(cmd))
        sys.exit(_stream(cmd, cwd=project_dir))

    # Render the SLURM profile from the pipeline config.
    profile_dir = project_dir / ".convgeno" / "slurm_profile"
    profile_dir.mkdir(parents=True, exist_ok=True)
    profile = render_slurm_profile(config, jobs)
    with open(profile_dir / "config.yaml", "w", encoding="utf-8") as f:
        yaml.dump(profile, f, default_flow_style=False, sort_keys=False)

    snakemake_cmd = _snakemake_common(snakefile, config_abs, mode) + [
        "--workflow-profile",
        str(profile_dir),
        "--jobs",
        str(jobs),
        "all",
    ]

    if local:
        print("Running Snakemake in the foreground (--local).")
        print("  " + " ".join(snakemake_cmd))
        sys.exit(_stream(snakemake_cmd, cwd=project_dir))

    # Default: submit a small orchestrator job that runs Snakemake.
    if config.runtime is None:
        print(
            "ERROR: No runtime configuration in the config. Run 'convgeno init' "
            "to detect conda paths (needed to activate the env in the "
            "orchestrator job)."
        )
        sys.exit(1)

    (project_dir / "logs").mkdir(parents=True, exist_ok=True)
    script_text = _orchestrator_script(config, snakemake_cmd, project_dir)
    script_dir = project_dir / ".convgeno"
    script_dir.mkdir(parents=True, exist_ok=True)
    script_path = script_dir / "convgeno_run.sh"
    script_path.write_text(script_text, encoding="utf-8")
    script_path.chmod(script_path.stat().st_mode | 0o755)

    if not skip_confirm:
        print(f"\n--- Orchestrator job script: {script_path} ---")
        print(script_text)
        print("--- End of script ---\n")
        response = input("Submit this orchestrator job? [Y/n]: ").strip().lower()
        if response.startswith("n"):
            print(f"Aborted. Script saved to: {script_path}")
            return

    # Reuse the OrthoFinder submitter for its automatic account fallback.
    from convgeno.cli.orthofinder_cmd import submit_sbatch

    try:
        sr = submit_sbatch(script_path)
    except FileNotFoundError:
        print(
            "ERROR: sbatch not found. Are you on a SLURM login node? "
            "Use 'convgeno run --local' inside an interactive allocation instead."
        )
        sys.exit(1)
    except RuntimeError as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("ERROR: sbatch submission timed out.")
        sys.exit(1)

    print("Orchestrator job submitted — it will run the whole pipeline as a DAG.")
    print(f"  Job ID:  {sr.job_id}")
    print(f"  Logs:    {project_dir}/logs/convgeno_run_{sr.job_id}.out")
    print("  Monitor: squeue -u $USER")
    print(f"  Cancel:  scancel {sr.job_id}")


def _stream(cmd: list[str], cwd: Path) -> int:
    """Run *cmd* inheriting stdio; return its exit code."""
    try:
        return subprocess.run(cmd, cwd=str(cwd)).returncode
    except FileNotFoundError:
        print(
            "ERROR: 'snakemake' not found. Install the environment with "
            "scripts/setup_env.sh and activate it (conda activate convgeno)."
        )
        return 1
