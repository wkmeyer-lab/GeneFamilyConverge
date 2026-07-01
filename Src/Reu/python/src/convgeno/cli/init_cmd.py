"""Interactive ``convgeno init`` wizard for first-time cluster setup."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig, normalize_optional_account
from convgeno.slurm.discovery import (
    detect_node_cpus,
    detect_partition_memory,
    detect_scratch_dir,
    discover_partitions,
    recommend_memory_mb,
)
from convgeno.slurm.runtime import (
    CondaRuntimeConfig,
    detect_conda_runtime,
    runtime_config_to_dict,
)


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


def _prompt_optional_int(message: str, default: int | None = None) -> int | None:
    """Prompt for an optional integer value.

    Returns *default* if the user presses Enter, ``None`` if they type
    nothing and *default* is ``None``.
    """
    if default is not None:
        prompt_str = f"{message} [{default}]: "
    else:
        prompt_str = f"{message} (press Enter to skip): "
    raw = input(prompt_str).strip()
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        print("  Please enter a valid integer.")
        return _prompt_optional_int(message, default)


def _derive_orthofinder_threads(recommended_physical: int) -> tuple[int, int]:
    """Derive OrthoFinder search and analysis thread counts."""
    if recommended_physical <= 0:
        return 16, 4
    return recommended_physical, max(recommended_physical // 4, 1)


def _parse_aligner_choice(user_input: str) -> str:
    """Parse an interactive OrthoFinder MSA aligner choice."""
    choice = user_input.strip().lower()
    if choice in {"", "1", "mafft"}:
        return "mafft"
    if choice in {"2", "famsa"}:
        return "famsa"
    return "mafft"


def _default_orthofinder_output_dir(
    project_dir: Path,
    *,
    timestamp: str | None = None,
) -> Path:
    """Return a fresh timestamped default output path for single-node OrthoFinder."""
    run_timestamp = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    return (
        project_dir
        / "Data"
        / "processed"
        / f"orthofinder_single_{run_timestamp}"
    )


def _format_detected_memory(value: int | None) -> str:
    return f"{value} MB" if value is not None else "unavailable"


def _memory_recommendation_basis(
    cpus_per_task: int,
    memory_detection: dict[str, int | None],
) -> str:
    max_mem_per_cpu = memory_detection["max_mem_per_cpu_mb"]
    def_mem_per_cpu = memory_detection["def_mem_per_cpu_mb"]
    min_node_memory = memory_detection["min_node_memory_mb"]
    if max_mem_per_cpu is not None and max_mem_per_cpu > 0:
        return f"MaxMemPerCPU={max_mem_per_cpu} MB x cpus_per_task={cpus_per_task}"
    if def_mem_per_cpu is not None and def_mem_per_cpu > 0:
        return f"DefMemPerCPU={def_mem_per_cpu} MB x cpus_per_task={cpus_per_task}"
    if min_node_memory is not None and min_node_memory > 0:
        return f"90% of minimum node memory ({min_node_memory} MB)"
    return "user-provided explicit memory request"


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

    cpu_detection = detect_node_cpus(partition_name)
    recommended_cpus = cpu_detection["recommended_physical"] or 16
    print(
        f"Detected {cpu_detection['node_count']} nodes in '{partition_name}', "
        f"{cpu_detection['min_cpus_per_node']} CPUs/node "
        f"(physical: {cpu_detection['physical_cores']}, "
        f"threads/core: {cpu_detection['threads_per_core']})"
    )
    print(
        f"Recommended CPUs per task: {recommended_cpus} "
        "(reserves 4 cores for memory headroom)"
    )

    cpus_per_task = _prompt_int("CPUs per task", default=recommended_cpus)
    memory_detection = detect_partition_memory(partition_name)
    print("Detected partition memory limits:")
    print(
        "  MaxMemPerCPU: "
        f"{_format_detected_memory(memory_detection['max_mem_per_cpu_mb'])}"
    )
    print(
        "  DefMemPerCPU: "
        f"{_format_detected_memory(memory_detection['def_mem_per_cpu_mb'])}"
    )
    print(
        "  Minimum node memory: "
        f"{_format_detected_memory(memory_detection['min_node_memory_mb'])}"
    )
    memory_basis = _memory_recommendation_basis(cpus_per_task, memory_detection)
    try:
        recommended_memory_mb = recommend_memory_mb(
            cpus_per_task=cpus_per_task,
            max_mem_per_cpu_mb=memory_detection["max_mem_per_cpu_mb"],
            def_mem_per_cpu_mb=memory_detection["def_mem_per_cpu_mb"],
            min_node_memory_mb=memory_detection["min_node_memory_mb"],
        )
        recommended_memory = f"{recommended_memory_mb}M"
        print(f"Recommended memory request: {recommended_memory}")
        print(f"Reason: {memory_basis}")
        memory_request = _prompt("Memory request", default=recommended_memory)
        if memory_request != recommended_memory:
            memory_basis = "user-provided explicit memory request"
    except ValueError:
        print(
            "Could not detect partition memory limits. Please enter an "
            "explicit SLURM memory request."
        )
        memory_request = _prompt("Memory request (e.g. 350400M)")
        memory_basis = "user-provided explicit memory request"

    search_threads, analysis_threads = _derive_orthofinder_threads(cpus_per_task)
    print(
        f"OrthoFinder threads: -t {search_threads} (sequence search), "
        f"-a {analysis_threads} (analysis)"
    )
    print("MSA aligner for OrthoFinder gene tree inference:")
    print(
        "  [1] mafft  — slower, more accurate "
        "(recommended: species tree used for state reconstruction)"
    )
    print(
        "  [2] famsa  — faster, slightly less accurate "
        "(use for large datasets where speed matters)"
    )
    msa_program = _parse_aligner_choice(_prompt("Select MSA aligner", default="1"))
    if msa_program == "famsa":
        print(
            "Using FAMSA. Ensure it is installed in your conda environment: "
            "conda install -c bioconda famsa"
        )
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

    # ## NEW: Detect scratch space for OrthoFinder's intermediate-file workload.
    scratch_detection = detect_scratch_dir()
    scratch_base = scratch_detection["scratch_base"]
    is_ephemeral_scratch = bool(scratch_detection["is_ephemeral"])
    if scratch_base is not None:
        scratch_kind = (
            "ephemeral — results will be copied back before job ends"
            if is_ephemeral_scratch
            else "persistent"
        )
        print(f"Detected scratch space: {scratch_base} ({scratch_kind})")
        if scratch_detection.get("write_granted"):
            print(
                f"  Added owner-write permission (chmod u+w) to {scratch_base} "
                "so it can be used as scratch."
            )
        print(
            "Using scratch will improve I/O performance for OrthoFinder's "
            "many intermediate files."
        )
    else:
        print(
            "No scratch space detected. OrthoFinder will run directly in the "
            "output directory."
        )
        print(
            "This may be slower on shared filesystems (e.g., Ceph) due to "
            "I/O pressure from intermediate files."
        )

    time_limit = _prompt("Job time limit (HH:MM:SS or D-HH:MM:SS)", default="72:00:00")
    if ":" not in time_limit:
        print(
            "  Warning: time limit format may be invalid. "
            "Expected HH:MM:SS or D-HH:MM:SS."
        )

    mail_user = _prompt_optional("Email for SLURM job notifications")
    account = normalize_optional_account(
        _prompt_optional("SLURM allocation/project account [optional, press Enter to omit]")
    )
    open_file_limit = _prompt_optional_int(
        "Maximum open files per job [optional, default 8192; press Enter to use default]",
        default=8192,
    )
    print(
        "  Lower analysis threads reduce open-file/shared-memory pressure "
        "during large OrthoFinder resume jobs."
    )
    print(
        "  Set orthofinder.analysis_threads in the config to override the "
        "auto-selected value."
    )

    slurm = SlurmConfig(
        partition=partition_name,
        time_limit=time_limit,
        cpus_per_task=cpus_per_task,
        mem=memory_request,
        mem_per_cpu=None,
        mail_user=mail_user,
        account=account,
        open_file_limit=open_file_limit,
        scratch_dir=str(scratch_base) if scratch_base is not None else None,
        is_ephemeral_scratch=is_ephemeral_scratch,
    )
    project_path = Path(project_dir)
    orthofinder = OrthoFinderConfig(
        input_dir=str(project_path / "Data/interim/cleaned_proteomes"),
        output_dir=str(_default_orthofinder_output_dir(project_path)),
        search_threads=search_threads,
        analysis_threads=analysis_threads,
        msa_program=msa_program,
    )

    # ---- Detect conda runtime configuration ----
    print("\n=== Detecting conda runtime ===\n")
    try:
        detected_runtime = detect_conda_runtime()
    except RuntimeError as e:
        print(f"WARNING: Could not detect conda runtime: {e}")
        print("You can configure runtime paths manually in the generated config file.")
        detected_runtime = None

    runtime_config: CondaRuntimeConfig | None = None
    if detected_runtime is not None:
        print("Detected conda runtime:")
        print(f"  Module: {detected_runtime.conda_module or '(none)'}")
        print(f"  Conda base: {detected_runtime.conda_base}")
        print(f"  Environment: {detected_runtime.conda_env_prefix}")
        print("")
        use_detected = _prompt("Use these settings? [Y/n]", default="Y")
        if use_detected.lower().startswith("n"):
            conda_module_str = _prompt_optional(
                "Conda module name (e.g. miniforge3/24.3.0-0)"
            )
            conda_base_str = _prompt(
                "Absolute path to conda installation",
                default=str(detected_runtime.conda_base),
            )
            conda_env_str = _prompt(
                "Absolute path to conda environment",
                default=str(detected_runtime.conda_env_prefix),
            )
            runtime_config = CondaRuntimeConfig(
                conda_module=conda_module_str,
                conda_base=Path(conda_base_str),
                conda_env_prefix=Path(conda_env_str),
            )
        else:
            runtime_config = detected_runtime

    config = PipelineConfig(
        project_dir=project_dir,
        conda_env=conda_env,
        slurm=slurm,
        orthofinder=orthofinder,
        runtime=runtime_config,
    )

    config.save(path)

    print(f"\nConfig written to {path}\n")
    print("Summary:")
    print(f"  Project directory:  {project_dir}")
    print(f"  Conda environment:  {conda_env}")
    print(f"  SLURM partition:    {partition_name}")
    print(f"  CPUs per task:      {cpus_per_task}")
    print(f"  Time limit:         {time_limit}")
    print(f"  Memory:             {memory_request}")
    print(f"  Memory basis:       {memory_basis}")
    if mail_user is not None:
        print(f"  Mail user:          {mail_user}")
    print(f"  SLURM account:      {account if account is not None else 'omitted'}")
    print(f"  Open-file limit:    {open_file_limit if open_file_limit is not None else 'not set'}")
    print(f"  Scratch directory:  {scratch_base if scratch_base is not None else 'not set'}")
    print(f"  OrthoFinder -t:     {search_threads}")
    print(f"  OrthoFinder -a:     {analysis_threads}")
    print(f"  OrthoFinder MSA:    {msa_program}")
    if runtime_config is not None:
        print(f"  Conda module:       {runtime_config.conda_module or '(none)'}")
        print(f"  Conda base:         {runtime_config.conda_base}")
        print(f"  Conda env prefix:   {runtime_config.conda_env_prefix}")
