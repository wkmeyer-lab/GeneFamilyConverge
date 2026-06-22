"""Interactive ``convgeno init`` wizard for first-time cluster setup."""

from __future__ import annotations

from pathlib import Path

from convgeno.slurm.config import PipelineConfig, SlurmConfig, normalize_optional_account
from convgeno.slurm.discovery import discover_partitions


def _prompt(message: str, default: str | None = None) -> str:
    """Prompt the user for input with an optional default value."""
    if default is not None:
        prompt_str = f"{message} [{default}]: "
    else:
        prompt_str = f"{message}: "

    value = input(prompt_str).strip()

    if not value and default is not None:
        return default
    if not value:
        print("  This field is required.")
        return _prompt(message, default)
    return value


def _prompt_optional(message: str) -> str | None:
    """Prompt for an optional value. Returns None if the user presses Enter."""
    prompt_str = f"{message} (press Enter to skip): "
    value = input(prompt_str).strip()
    return value if value else None


def _prompt_int(message: str, default: int) -> int:
    """Prompt for an integer value with a default."""
    raw = _prompt(message, default=str(default))
    try:
        value = int(raw)
    except ValueError:
        print("  Please enter a valid integer.")
        return _prompt_int(message, default)
    if value < 1:
        print("  Value must be at least 1.")
        return _prompt_int(message, default)
    return value


def run_init(output_path: str = "pipeline_config.yaml") -> None:
    """Run the interactive init wizard to create pipeline_config.yaml."""
    path = Path(output_path)

    if path.exists():
        answer = input(
            f"Config file {output_path} already exists. Overwrite? [y/N]: "
        ).strip()
        if not answer.lower().startswith("y"):
            print("Aborted.")
            return

    print("\n=== convgeno init ===\n")
    print("This will create your pipeline configuration file.\n")

    project_dir = _prompt(
        "Project directory on the cluster filesystem", default=str(Path.cwd())
    )
    conda_env = _prompt("Conda environment name", default="convgeno")

    partitions = discover_partitions()
    selected_partition = None

    if partitions:
        print("\nDetected SLURM partitions:")
        for i, part in enumerate(partitions):
            print(f"  [{i + 1}] {part}")

        default_index = 1
        for i, part in enumerate(partitions):
            if part.is_default:
                default_index = i + 1
                break

        while True:
            choice_str = _prompt("Select partition number", default=str(default_index))
            try:
                choice = int(choice_str)
                if 1 <= choice <= len(partitions):
                    selected_partition = partitions[choice - 1]
                    partition_name = selected_partition.name
                    break
                print(f"  Please enter a number between 1 and {len(partitions)}.")
            except ValueError:
                print("  Please enter a valid number.")
    else:
        print("\nNo SLURM partitions detected (sinfo not available).")
        partition_name = _prompt("Enter your SLURM partition name")

    cpus_per_task = _prompt_int("CPUs per task", default=16)
    if (
        selected_partition is not None
        and selected_partition.max_cpus_per_node > 0
        and cpus_per_task > selected_partition.max_cpus_per_node
    ):
        print(
            f"  Warning: requested {cpus_per_task} CPUs but partition "
            f"'{partition_name}' has {selected_partition.max_cpus_per_node} "
            f"CPUs per node."
        )

    time_limit = _prompt("Job time limit (HH:MM:SS or D-HH:MM:SS)", default="48:00:00")
    if ":" not in time_limit:
        print(
            "  Warning: time limit format may be invalid. "
            "Expected HH:MM:SS or D-HH:MM:SS."
        )

    mem_per_cpu = _prompt("Memory per CPU", default="4G")
    mail_user = _prompt_optional("Email for SLURM job notifications")
    account = normalize_optional_account(
        _prompt_optional("SLURM allocation/project account [optional, press Enter to omit]")
    )

    slurm = SlurmConfig(
        partition=partition_name,
        time_limit=time_limit,
        cpus_per_task=cpus_per_task,
        mem_per_cpu=mem_per_cpu,
        mail_user=mail_user,
        account=account,
    )
    config = PipelineConfig(
        project_dir=project_dir,
        conda_env=conda_env,
        slurm=slurm,
    )

    config.save(path)

    print(f"\nConfig written to {path}\n")
    print("Summary:")
    print(f"  Project directory:  {project_dir}")
    print(f"  Conda environment:  {conda_env}")
    print(f"  SLURM partition:    {partition_name}")
    print(f"  CPUs per task:      {cpus_per_task}")
    print(f"  Time limit:         {time_limit}")
    print(f"  Memory per CPU:     {mem_per_cpu}")
    if mail_user is not None:
        print(f"  Mail user:          {mail_user}")
    print(f"  SLURM account:      {account if account is not None else 'omitted'}")
