"""
Tests for convgeno.slurm.script_generator

Run with:  pytest tests/test_script_generator.py -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig
from convgeno.slurm.runtime import CondaRuntimeConfig
from convgeno.slurm.script_generator import generate_orthofinder_script, write_script


@pytest.fixture()
def sample_runtime() -> CondaRuntimeConfig:
    return CondaRuntimeConfig(
        conda_module="miniforge3/24.3.0-0",
        conda_base=Path("/share/apps/miniforge3/24.3.0-0"),
        conda_env_prefix=Path("/home/prm526/.conda/envs/convgeno"),
    )


@pytest.fixture()
def sample_config(sample_runtime) -> PipelineConfig:
    return PipelineConfig(
        project_dir="/share/ceph/project",
        conda_env="convgeno",
        slurm=SlurmConfig(
            partition="hawkcpu",
            cpus_per_task=16,
            time_limit="72:00:00",
            mem="350400M",
        ),
        orthofinder=OrthoFinderConfig(
            input_dir="/share/ceph/project/Data/interim/cleaned_proteomes",
            output_dir="/share/ceph/project/Data/processed/orthofinder",
        ),
        runtime=sample_runtime,
    )


class TestGenerateOrthoFinderScript:
    def test_starts_with_shebang(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert script.startswith("#!/bin/bash\n")

    def test_contains_sbatch_directives(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "#SBATCH --partition=hawkcpu" in script
        assert "#SBATCH --time=72:00:00" in script
        assert "#SBATCH --cpus-per-task=16" in script
        assert "#SBATCH --mem=350400M" in script
        assert "#SBATCH --job-name=convgeno_orthofinder" in script

    def test_does_not_emit_mem_zero_unless_explicit(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=sample_runtime,
        )

        script = generate_orthofinder_script(config)

        assert "#SBATCH --mem=0" not in script

    def test_explicit_mem_zero_is_preserved(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu", mem="0"),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=sample_runtime,
        )

        script = generate_orthofinder_script(config)

        assert "#SBATCH --mem=0" in script

    def test_contains_conda_activation(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert 'conda activate "$CONDA_ENV"' in script

    def test_contains_orthofinder_command(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "orthofinder" in script
        assert 'INPUT_DIR="/share/ceph/project/Data/interim/cleaned_proteomes"' in script
        assert 'OUTPUT_DIR="/share/ceph/project/Data/processed/orthofinder"' in script
        assert '-f "$EFFECTIVE_INPUT"' in script
        assert '-o "$EFFECTIVE_OUTPUT"' in script
        assert 'ORTHOFINDER_SEARCH_THREADS="16"' in script
        assert 'ORTHOFINDER_ANALYSIS_THREADS="8"' in script
        assert '-t "$ORTHOFINDER_SEARCH_THREADS"' in script
        assert '-a "$ORTHOFINDER_ANALYSIS_THREADS"' in script
        assert "-S diamond" in script
        # OrthoFinder's -M takes a method, not a program name.
        assert "-M msa" in script
        assert 'ORTHOFINDER_MSA_PROGRAM="mafft"' in script
        assert '-A "$ORTHOFINDER_MSA_PROGRAM"' in script
        assert "-T fasttree" in script
        assert "-M mafft" not in script

    def test_contains_input_validation(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert 'if [ ! -d "$INPUT_DIR" ]' in script
        assert "FASTA_COUNT" in script

    def test_does_not_create_output_dir_directly(self, sample_config: PipelineConfig):
        # OrthoFinder refuses to run with an existing -o output directory.
        # The script must create only the PARENT directory.
        script = generate_orthofinder_script(sample_config)
        assert 'mkdir -p "$(dirname "$OUTPUT_DIR")"' in script
        assert 'mkdir -p "$OUTPUT_DIR"' not in script

    def test_fails_loudly_if_output_dir_exists(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert 'if [ -e "$OUTPUT_DIR" ]' in script
        assert "already exists" in script

    def test_uses_help_for_version_probe_not_unsupported_flag(
        self, sample_config: PipelineConfig
    ):
        # OrthoFinder 2.5.5 does not support --version. Use -h instead.
        script = generate_orthofinder_script(sample_config)
        assert "orthofinder --version" not in script
        assert "orthofinder -h" in script

    def test_contains_exit_code_check(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "EXIT_CODE=$?" in script
        assert 'exit "$EXIT_CODE"' in script

    def test_disables_errexit_only_for_orthofinder_command(
        self, sample_config: PipelineConfig
    ):
        script = generate_orthofinder_script(sample_config)

        set_plus_index = script.index("set +e")
        command_index = script.index("orthofinder \\\n", set_plus_index)
        exit_code_index = script.index("EXIT_CODE=$?", command_index)
        set_minus_index = script.index("set -e", exit_code_index)

        assert set_plus_index < command_index < exit_code_index < set_minus_index

    def test_failure_handler_exits_with_orthofinder_status(
        self, sample_config: PipelineConfig
    ):
        script = generate_orthofinder_script(sample_config)

        assert 'if [ "$EXIT_CODE" -ne 0 ]; then' in script
        assert "ERROR: OrthoFinder exited with code $EXIT_CODE" in script
        assert 'exit "$EXIT_CODE"' in script

    def test_contains_timing(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "START_SECONDS" in script
        assert "Duration:" in script

    def test_set_euo_pipefail(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)
        assert "set -euo pipefail" in script

    def test_raises_without_orthofinder_config(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            runtime=sample_runtime,
        )
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_orthofinder_script(config)

    def test_with_extra_orthofinder_args(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(
                input_dir="/in",
                output_dir="/out",
                extra_args=["--fewer-files"],
            ),
            runtime=sample_runtime,
        )
        script = generate_orthofinder_script(config)
        assert "--fewer-files" in script

    def test_generated_script_contains_export_all(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        assert "#SBATCH --export=ALL" in script

    def test_generated_script_no_conda_info_base_as_primary(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        primary_end = script.index("elif command -v conda")
        assert "conda info --base" not in script[:primary_end]

    def test_contains_bootstrap_markers(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        assert "# ---- convgeno runtime bootstrap ----" in script
        assert "# ---- end convgeno runtime bootstrap ----" in script

    def test_sources_conda_sh_with_absolute_path(self, sample_config):
        script = generate_orthofinder_script(sample_config)
        assert '[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]' in script
        assert 'source "$CONDA_BASE/etc/profile.d/conda.sh"' in script

    def test_raises_without_runtime(self):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=None,
        )
        with pytest.raises(ValueError, match="runtime"):
            generate_orthofinder_script(config)

    def test_explicit_runtime_overrides_config(self):
        alt_runtime = CondaRuntimeConfig(
            conda_module=None,
            conda_base=Path("/alt/conda"),
            conda_env_prefix=Path("/alt/envs/myenv"),
        )
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=CondaRuntimeConfig(
                conda_module="old/1.0",
                conda_base=Path("/old/conda"),
                conda_env_prefix=Path("/old/envs/convgeno"),
            ),
        )
        script = generate_orthofinder_script(config, alt_runtime)
        assert 'CONDA_BASE="/alt/conda"' in script
        assert 'conda activate "$CONDA_ENV"' in script

    def test_script_without_scratch(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)

        assert "JOB_SCRATCH" not in script
        assert "rsync" not in script
        assert "SCRATCH_INPUT" not in script
        assert "cleanup_scratch" not in script

    def test_script_with_scratch(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                scratch_dir="/share/ceph/scratch/testuser",
            ),
            orthofinder=OrthoFinderConfig(
                input_dir="/share/ceph/project/Data/interim/cleaned_proteomes",
                output_dir="/share/ceph/project/Data/processed/orthofinder",
            ),
            runtime=sample_runtime,
        )

        script = generate_orthofinder_script(config)

        # The preferred base is baked in, but JOB_SCRATCH derives from the
        # runtime-resolved base so the script can self-heal / fall back.
        assert 'PREFERRED_SCRATCH="/share/ceph/scratch/testuser"' in script
        assert 'JOB_SCRATCH="$SCRATCH_BASE/${SLURM_JOB_ID}"' in script
        assert "SCRATCH_INPUT=" in script
        assert "SCRATCH_OUTPUT=" in script
        assert "export TMPDIR=" in script
        assert script.count("rsync -a") >= 2
        assert "cleanup_scratch" in script
        assert "trap cleanup_scratch" in script
        assert 'orthofinder \\\n  -f "$EFFECTIVE_INPUT" \\\n  -o "$EFFECTIVE_OUTPUT"' in script
        assert "orthofinder -f /share/ceph/project/Data/interim/cleaned_proteomes" not in script

    def test_scratch_runtime_resolver_present(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                scratch_dir="/share/ceph/scratch/testuser",
            ),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=sample_runtime,
        )

        script = generate_orthofinder_script(config)

        # Self-healing resolver: chmod u+w when owned-but-unwritable, real
        # write probe, then fall back to /tmp/scratch.
        assert "resolve_scratch_base()" in script
        assert "chmod u+w" in script
        assert ".convgeno_wtest" in script
        assert 'FALLBACK_SCRATCH="/tmp/scratch"' in script
        assert 'SCRATCH_BASE="$(resolve_scratch_base "$PREFERRED_SCRATCH")"' in script

    def test_scratch_resolver_absent_without_scratch(
        self, sample_config: PipelineConfig
    ):
        script = generate_orthofinder_script(sample_config)

        assert "resolve_scratch_base" not in script
        assert "PREFERRED_SCRATCH" not in script
        assert "chmod u+w" not in script

    def test_scratch_cleanup_present(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                scratch_dir="/share/ceph/scratch/testuser",
            ),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out"),
            runtime=sample_runtime,
        )

        script = generate_orthofinder_script(config)

        assert 'rm -rf "$JOB_SCRATCH"' in script

    def test_effective_paths_without_scratch(self, sample_config: PipelineConfig):
        script = generate_orthofinder_script(sample_config)

        assert 'EFFECTIVE_INPUT="$INPUT_DIR"' in script
        assert 'EFFECTIVE_OUTPUT="$OUTPUT_DIR"' in script
        assert '-f "$EFFECTIVE_INPUT"' in script
        assert '-o "$EFFECTIVE_OUTPUT"' in script


class TestWriteScript:
    def test_creates_file(self, tmp_path: Path):
        script_path = tmp_path / "test_job.sh"
        write_script("#!/bin/bash\necho hello", script_path)
        assert script_path.exists()
        assert script_path.read_text(encoding="utf-8") == "#!/bin/bash\necho hello"

    @pytest.mark.skipif(
        sys.platform == "win32",
        reason="Windows does not support Unix file permission bits",
    )
    def test_is_executable(self, tmp_path: Path):
        script_path = tmp_path / "job.sh"
        write_script("#!/bin/bash\necho test", script_path)
        assert script_path.stat().st_mode & 0o111

    def test_creates_parent_dirs(self, tmp_path: Path):
        script_path = tmp_path / "subdir" / "nested" / "job.sh"
        write_script("#!/bin/bash", script_path)
        assert script_path.exists()
