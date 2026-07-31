"""Interactive ``convgeno init`` wizard for first-time cluster setup."""

from __future__ import annotations

import dataclasses
import re
from datetime import datetime
from pathlib import Path

from convgeno.external.config import OrthoFinderConfig
from convgeno.io.fasta import FASTA_EXTENSIONS
from convgeno.slurm.config import (
    PhenotypeTreeConfig,
    PipelineConfig,
    SlurmConfig,
    SpeciesTreeConfig,
    UltrametricConfig,
    normalize_optional_account,
)
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
from convgeno.validation.trees import (
    check_tips_match_species,
    is_ultrametric,
    root_to_tip_depths,
    ultrametric_deviation,
    validate_tree,
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


def _prompt_float(message: str) -> float:
    """Prompt for a positive float value (re-prompts until valid)."""
    while True:
        raw = _prompt(message)
        try:
            value = float(raw)
        except ValueError:
            print("  Please enter a number (e.g. 94 or 94.0).")
            continue
        if value <= 0:
            print("  Value must be greater than 0.")
            continue
        return value


def _prompt_yes_no(message: str, default: bool = False) -> bool:
    """Prompt for a yes/no answer, returning *default* on empty input."""
    suffix = "[Y/n]" if default else "[y/N]"
    raw = input(f"{message} {suffix}: ").strip().lower()
    if not raw:
        return default
    return raw.startswith("y")


def _prompt_int_optional(message: str) -> int | None:
    """Prompt for an optional positive integer. Returns None on empty input."""
    while True:
        raw = input(f"{message} (press Enter to skip): ").strip()
        if not raw:
            return None
        try:
            value = int(raw)
        except ValueError:
            print("  Please enter a valid integer.")
            continue
        if value < 1:
            print("  Value must be at least 1.")
            continue
        return value


#: A species name is used verbatim as an OrthoFinder tip label and embedded in
#: the r8s ``--calibration NAME:SP1,SP2:AGE`` argument, so it must not contain
#: the ``:``/``,`` delimiters or shell-hostile characters.
_SAFE_SPECIES_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def _proteome_species_names(input_dir: Path) -> set[str]:
    """Return the species names (FASTA basenames) present in *input_dir*.

    Empty if the directory does not exist yet (proteomes not prepared), in which
    case calibration species can't be validated at ``init`` time.
    """
    if not input_dir.is_dir():
        return set()
    return {
        entry.stem
        for entry in input_dir.iterdir()
        if entry.is_file() and entry.suffix.lower() in FASTA_EXTENSIONS
    }


def _is_valid_species(name: str, available: set[str]) -> bool:
    """True if *name* is well-formed and (when known) present in *available*."""
    if not _SAFE_SPECIES_RE.match(name):
        return False
    return not available or name in available


def _report_bad_species(name: str, available: set[str]) -> None:
    if not _SAFE_SPECIES_RE.match(name):
        print(
            f"  '{name}' has unsupported characters; use only letters, digits, "
            "'.', '_', '-' (the proteome basename)."
        )
    elif available:
        sample = ", ".join(sorted(available)[:8])
        print(f"  '{name}' is not one of your proteome species. Examples: {sample}")


def _prompt_species(message: str, available: set[str]) -> str:
    """Prompt (required) for a valid species name, re-prompting on error."""
    while True:
        name = _prompt(message)
        if _is_valid_species(name, available):
            return name
        _report_bad_species(name, available)


def _prompt_calibration(input_dir: Path) -> UltrametricConfig:
    """Prompt for the two-species + divergence-time r8s calibration.

    Species are validated against the proteome files in *input_dir* when those
    exist. Pressing Enter at the first prompt skips calibration (the tree is
    then made ultrametric in relative time).
    """
    print("\n=== Species tree calibration (r8s ultrametric step) ===")
    print(
        "r8s scales the OrthoFinder species tree to time from ONE calibration:\n"
        "two species and their divergence time in millions of years. Name each\n"
        "species EXACTLY as its proteome file basename (no extension), e.g.\n"
        "'Homo_sapiens'."
    )
    available = _proteome_species_names(input_dir)
    if available:
        print(f"  {len(available)} proteome species found in {input_dir}.")
    else:
        print(
            f"  (No proteomes in {input_dir} yet — names will be verified against "
            "the species tree at run time.)"
        )

    first = _prompt_optional(
        "First calibration species (Enter to skip and use relative time)"
    )
    if first is None:
        print(
            "  No calibration set: the tree will be made ultrametric in RELATIVE "
            "time (root-anchored). Re-run 'convgeno init' to add one later."
        )
        return UltrametricConfig()

    if not _is_valid_species(first, available):
        _report_bad_species(first, available)
        first = _prompt_species("First calibration species", available)

    second = _prompt_species("Second calibration species", available)
    while second == first:
        print("  The two species must be different.")
        second = _prompt_species("Second calibration species", available)

    age = _prompt_float("Divergence time between them (millions of years)")
    print(f"  Calibration set: MRCA({first}, {second}) = {age:g} Myr.")
    return UltrametricConfig(species_a=first, species_b=second, divergence_my=age)


def _verify_ultrametric(tree_path: Path) -> bool:
    """Verify a tree is ultrametric, reporting the numeric root-to-tip spread.

    Returns True only when :func:`is_ultrametric` confirms it. On failure the
    per-tip root-to-tip depths are printed so the user can see why, and the tree
    is treated as non-ultrametric (r8s will date it), since CAFE-5 needs a
    genuinely ultrametric tree.
    """
    depths = root_to_tip_depths(tree_path)
    deviation = ultrametric_deviation(tree_path)
    height = max(depths.values()) if depths else 0.0
    relative = (deviation / height) if height > 0 else 0.0
    print(
        f"  Root-to-tip depth: height={height:g}, max-min deviation="
        f"{deviation:g} ({relative:.2%} of height)."
    )
    if is_ultrametric(tree_path):
        print("  Verified ultrametric: the r8s dating step will be skipped.")
        return True
    print(
        "  This tree is NOT ultrametric within tolerance. It will be treated as\n"
        "  non-ultrametric so r8s can date it (CAFE-5 requires an ultrametric tree)."
    )
    print("  Per-tip root-to-tip depths:")
    for name, depth in sorted(depths.items(), key=lambda kv: kv[1]):
        print(f"    {name}: {depth:g}")
    return False


def _prompt_species_tree(input_dir: Path) -> SpeciesTreeConfig | None:
    """Prompt for an optional user-supplied species tree.

    Returns ``None`` when the user declines (OrthoFinder infers the tree).
    Otherwise the Newick file is validated (parseable, rooted, binary),
    ultrametricity is verified when the user claims it, and — when proteomes are
    already staged — the tip labels are checked against the proteome species for
    instant feedback. Re-prompts on an unreadable or invalid tree; pressing
    Enter at the re-prompt backs out to the OrthoFinder default.
    """
    print("\n=== Species tree (optional: bring your own) ===")
    print(
        "By default OrthoFinder infers the species tree. If you already have a\n"
        "trusted tree for these proteomes, provide it here as a Newick file. Tip\n"
        "labels must match your proteome filenames (no extension)."
    )
    while True:
        path_str = _prompt_optional(
            "Path to your own species tree (Newick) [Enter to use OrthoFinder's]"
        )
        if path_str is None:
            return None

        tree_path = Path(path_str).expanduser()
        if not tree_path.is_file():
            print(f"  File not found: {tree_path}")
            continue

        errors = validate_tree(tree_path)
        if errors:
            print("  This tree is not usable by the pipeline:")
            for err in errors:
                print(f"    - {err}")
            continue

        available = _proteome_species_names(input_dir)
        if available:
            in_tree, in_proteomes = check_tips_match_species(
                tree_path, sorted(available)
            )
            if in_proteomes:
                sample = ", ".join(sorted(in_proteomes)[:8])
                print(
                    f"  Warning: {len(in_proteomes)} proteome species missing from "
                    f"the tree: {sample}"
                )
            if in_tree:
                sample = ", ".join(sorted(in_tree)[:8])
                print(
                    f"  Warning: {len(in_tree)} tree tips are not among your "
                    f"proteomes: {sample}"
                )
            if in_tree or in_proteomes:
                if not _prompt_yes_no("Continue with this tree anyway?", default=True):
                    continue
            else:
                print(f"  Tip labels match all {len(available)} proteome species.")
        else:
            print(
                f"  (No proteomes in {input_dir} yet — tip labels will be checked "
                "against the species set at run time.)"
            )

        verified_ultra = False
        if _prompt_yes_no(
            "Is this tree ultrametric (time-calibrated)?", default=False
        ):
            verified_ultra = _verify_ultrametric(tree_path)

        num_sites = None
        if not verified_ultra:
            num_sites = _prompt_int_optional(
                "Number of alignment sites behind the tree's branch lengths "
                "(for r8s) [Enter to use OrthoFinder's alignment length]"
            )

        return SpeciesTreeConfig(
            path=str(tree_path.resolve()),
            is_ultrametric=verified_ultra,
            num_sites=num_sites,
        )


def _read_tsv_column(path: Path, column: str) -> set[str]:
    """Return the set of non-empty values in a named TSV column (best-effort)."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return set()
    if not lines:
        return set()
    header = lines[0].split("\t")
    if column not in header:
        return set()
    idx = header.index(column)
    values: set[str] = set()
    for line in lines[1:]:
        cells = line.split("\t")
        if idx < len(cells):
            value = cells[idx].strip()
            if value:
                values.add(value)
    return values


def _prompt_phenotype_tree(input_dir: Path) -> PhenotypeTreeConfig | None:
    """Prompt for an optional phenotype table for the automatic CAFE -y tree.

    Returns ``None`` when the user declines (no phenotype tree is built).
    Otherwise the TSV is checked to exist and to hold the id/phenotype columns,
    and — when proteomes are already staged — phenotype coverage of the species
    set is reported. Pressing Enter at the first prompt skips the feature.
    """
    print(
        "\n=== Phenotype tip data (optional: automatic categorical phenotype tree) ==="
    )
    print(
        "Provide a phenotype table and the pipeline will AUTOMATICALLY build the\n"
        "categorical phenotype tree (the CAFE -y multi-lambda tree) right after the\n"
        "species tree is made ultrametric. The table is a TSV with a column of tip\n"
        "labels (= proteome filenames, no extension) and a phenotype column."
    )
    while True:
        path_str = _prompt_optional(
            "Path to your phenotype table (TSV) [Enter to skip]"
        )
        if path_str is None:
            return None

        table_path = Path(path_str).expanduser()
        if not table_path.is_file():
            print(f"  File not found: {table_path}")
            continue

        try:
            first_line = table_path.read_text(encoding="utf-8").splitlines()[0]
        except (OSError, IndexError):
            print(f"  Could not read a header row from {table_path}")
            continue
        columns = first_line.rstrip("\n").split("\t")
        if len(columns) < 2:
            print(
                "  This does not look tab-separated (< 2 columns in the header). "
                "Provide a TSV."
            )
            continue

        id_col = _prompt("Tip-label column name", default="species")
        pheno_col = _prompt("Phenotype column name", default="phenotype")
        missing_cols = [c for c in (id_col, pheno_col) if c not in columns]
        if missing_cols:
            print(
                f"  Column(s) not in the table header: {', '.join(missing_cols)}\n"
                f"  Available columns: {', '.join(columns)}"
            )
            if not _prompt_yes_no("Use this table anyway?", default=False):
                continue

        available = _proteome_species_names(input_dir)
        if available and id_col in columns:
            labelled = _read_tsv_column(table_path, id_col)
            covered = available & labelled
            missing = available - labelled
            print(
                f"  {len(covered)}/{len(available)} proteome species have a "
                "phenotype in this table."
            )
            if missing:
                sample = ", ".join(sorted(missing)[:8])
                print(
                    f"  {len(missing)} without one (they will get a 'background' "
                    f"class): {sample}"
                )

        model = _prompt("ASR rate model — ER, SYM, or ARD", default="ER")
        model = model.strip().upper()
        if model not in {"ER", "SYM", "ARD"}:
            print(f"  Unrecognized model '{model}'; using ER.")
            model = "ER"

        return PhenotypeTreeConfig(
            table=str(table_path.resolve()),
            id_col=id_col,
            pheno_col=pheno_col,
            model=model,
        )


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
    mode: str,
    *,
    timestamp: str | None = None,
) -> Path:
    """Return a fresh timestamped OrthoFinder output path, named by ``mode``.

    The execution mode (``"multinode"`` / ``"singlenode"``) is baked into the
    directory name (``orthofinder_<mode>_<ts>``) so single-node and multi-node
    runs never collide on the same ``-o`` directory and their results are easy
    to tell apart.
    """
    run_timestamp = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    return (
        project_dir
        / "Data"
        / "processed"
        / f"orthofinder_{mode}_{run_timestamp}"
    )


def _mode_config_path(base: Path, mode: str) -> Path:
    """Mode-specific config filename derived from *base*.

    ``pipeline_config.yaml`` + ``multinode`` -> ``pipeline_config_multinode.yaml``.
    An existing ``_multinode`` / ``_singlenode`` suffix on *base* is stripped
    first so re-running ``init`` never doubles it.
    """
    stem = base.stem
    for suffix in ("_multinode", "_singlenode"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return base.with_name(f"{stem}_{mode}{base.suffix or '.yaml'}")


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
    """Run the interactive init wizard, writing one config per execution mode.

    Writes ``pipeline_config_multinode.yaml`` and ``pipeline_config_singlenode.yaml``
    (derived from *output_path*), each with a mode-named ``output_dir``, so the
    two jobs can be launched together without colliding on OrthoFinder's ``-o``.
    """
    base_path = Path(output_path)
    mode_paths = {
        m: _mode_config_path(base_path, m) for m in ("multinode", "singlenode")
    }

    existing = [p for p in mode_paths.values() if p.exists()]
    if existing:
        listed = ", ".join(str(p) for p in existing)
        answer = input(
            f"Config file(s) {listed} already exist. Overwrite? [y/N]: "
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
    slurm = SlurmConfig(
        partition=partition_name,
        time_limit=time_limit,
        cpus_per_task=cpus_per_task,
        mem=memory_request,
        mem_per_cpu=None,
        mail_user=mail_user,
        account=account,
        scratch_dir=str(scratch_base) if scratch_base is not None else None,
        is_ephemeral_scratch=is_ephemeral_scratch,
    )
    project_path = Path(project_dir)
    run_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    # output_dir is filled in per mode when the two configs are written below.
    orthofinder = OrthoFinderConfig(
        input_dir=str(project_path / "Data/interim/cleaned_proteomes"),
        output_dir="",
        search_threads=search_threads,
        analysis_threads=analysis_threads,
        msa_program=msa_program,
    )

    species_tree = _prompt_species_tree(Path(orthofinder.input_dir))
    if species_tree is not None and species_tree.is_ultrametric:
        print(
            "\nUltrametric species tree supplied; skipping the r8s calibration "
            "(r8s will not run)."
        )
        ultrametric = UltrametricConfig()
    else:
        ultrametric = _prompt_calibration(Path(orthofinder.input_dir))

    phenotype_tree = _prompt_phenotype_tree(Path(orthofinder.input_dir))

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
        ultrametric=ultrametric,
        species_tree=species_tree,
        phenotype_tree=phenotype_tree,
    )

    # One config per execution mode, each with a mode-named output_dir sharing
    # this run's timestamp, so single-node and multi-node jobs never collide on
    # OrthoFinder's -o directory.
    written: list[tuple[str, Path, str]] = []
    for mode, mode_path in mode_paths.items():
        out_dir = str(
            _default_orthofinder_output_dir(project_path, mode, timestamp=run_ts)
        )
        mode_of = dataclasses.replace(config.orthofinder, output_dir=out_dir)
        mode_cfg = dataclasses.replace(config, orthofinder=mode_of)
        mode_cfg.save(mode_path)
        written.append((mode, mode_path, out_dir))

    print("\nConfigs written (one per execution mode):")
    for mode, mode_path, out_dir in written:
        print(f"  {mode:10s} -> {mode_path}")
        print(f"               output_dir: {out_dir}")
    print("")
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
    print(f"  Scratch directory:  {scratch_base if scratch_base is not None else 'not set'}")
    print(f"  OrthoFinder -t:     {search_threads}")
    print(f"  OrthoFinder -a:     {analysis_threads}")
    print(f"  OrthoFinder MSA:    {msa_program}")
    if species_tree is not None and species_tree.has_tree():
        kind = (
            "ultrametric (r8s skipped)"
            if species_tree.is_ultrametric
            else "non-ultrametric (r8s will date it)"
        )
        print(f"  User species tree:  {species_tree.path}")
        print(f"                      {kind}")
        if species_tree.num_sites is not None:
            print(f"                      num_sites: {species_tree.num_sites}")
    else:
        print("  Species tree:       OrthoFinder-inferred")
    if ultrametric.has_calibration():
        print(
            f"  r8s calibration:    MRCA({ultrametric.species_a}, "
            f"{ultrametric.species_b}) = {ultrametric.divergence_my:g} Myr"
        )
    elif species_tree is not None and species_tree.is_ultrametric:
        print("  r8s calibration:    n/a (ultrametric tree supplied)")
    else:
        print("  r8s calibration:    none (relative-time, root-anchored)")
    if phenotype_tree is not None and phenotype_tree.has_table():
        print(f"  Phenotype tree:     auto CAFE -y from {phenotype_tree.table}")
        print(
            f"                      cols {phenotype_tree.id_col} -> "
            f"{phenotype_tree.pheno_col}, model {phenotype_tree.model}"
        )
    else:
        print("  Phenotype tree:     none (no phenotype table given)")
    if runtime_config is not None:
        print(f"  Conda module:       {runtime_config.conda_module or '(none)'}")
        print(f"  Conda base:         {runtime_config.conda_base}")
        print(f"  Conda env prefix:   {runtime_config.conda_env_prefix}")
