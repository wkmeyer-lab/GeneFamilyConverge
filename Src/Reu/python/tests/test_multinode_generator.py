"""
Tests for convgeno.slurm.multinode_generator

Run with:  pytest tests/test_multinode_generator.py -v
"""

from __future__ import annotations

import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from convgeno.external.config import OrthoFinderConfig
from convgeno.slurm.config import PipelineConfig, SlurmConfig
from convgeno.slurm.multinode_generator import (
    COMMAND_LINE_REGEX,
    COMMAND_LINE_REGEX_EXTRACT,
    _proteome_volume,
    derive_prepare_cpus,
    derive_prepare_memory_mb,
    derive_prepare_walltime,
    generate_prepare_script,
    generate_resume_script,
    generate_search_array_script,
)
from convgeno.slurm.runtime import CondaRuntimeConfig


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
            account="wym219",
        ),
        orthofinder=OrthoFinderConfig(
            input_dir="/share/ceph/project/proteomes",
            output_dir="/share/ceph/project/results",
            search_threads=16,
            analysis_threads=8,
        ),
        runtime=sample_runtime,
    )


class TestPrepareScript:
    def test_shebang(self, sample_config):
        assert generate_prepare_script(sample_config).startswith("#!/bin/bash\n")

    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_prepare" in generate_prepare_script(sample_config)

    def test_single_node(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert "--nodes=1" in script
        # CPU = min(16, cpus_per_task); sample_config has cpus_per_task=16.
        assert "--cpus-per-task=16" in script

    def test_short_time_limit(self, sample_config):
        assert "--time=02:00:00" in generate_prepare_script(sample_config)

    def test_runs_op_flag(self, sample_config):
        assert "-op" in generate_prepare_script(sample_config)

    def test_captures_commands(self, sample_config):
        assert "diamond_commands.txt" in generate_prepare_script(sample_config)

    def test_uses_partition_from_config(self, sample_config):
        assert "--partition=hawkcpu" in generate_prepare_script(sample_config)

    def test_includes_account(self, sample_config):
        assert "--account=wym219" in generate_prepare_script(sample_config)

    def test_does_not_create_output_dir_directly(self, sample_config):
        # OrthoFinder refuses to run with an existing -o directory.
        # Only the parent should be created.
        script = generate_prepare_script(sample_config)
        assert 'mkdir -p "$OUTPUT_DIR"' not in script
        assert 'mkdir -p "$OUTPUT_PARENT"' in script

    def test_fails_loudly_if_output_dir_exists(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert 'if [ -e "$OUTPUT_DIR" ]' in script
        assert "already exists" in script

    def test_aux_files_live_in_parent_not_in_output_dir(self, sample_config):
        # Bug: bash opens redirect targets before exec'ing the command,
        # so any "> $OUTPUT_DIR/foo" fails because OrthoFinder hasn't
        # created $OUTPUT_DIR yet. All aux files must live in
        # $OUTPUT_PARENT/${RUN_NAME}_*.
        script = generate_prepare_script(sample_config)
        assert (
            'PREPARE_LOG="$OUTPUT_PARENT/${RUN_NAME}_prepare_full_stdout.log"'
            in script
        )
        assert (
            'COMMANDS_FILE="$OUTPUT_PARENT/${RUN_NAME}_diamond_commands.txt"'
            in script
        )
        assert (
            'WORK_DIR_FILE="$OUTPUT_PARENT/${RUN_NAME}_working_dir_path.txt"'
            in script
        )
        assert '> "$PREPARE_LOG" 2>&1' in script

    def test_no_writes_into_output_dir_before_orthofinder(self, sample_config):
        # Belt and suspenders: assert none of the previously buggy
        # paths reappear via a copy/paste regression.
        script = generate_prepare_script(sample_config)
        assert '"$OUTPUT_DIR/prepare_full_stdout.log"' not in script
        assert '"$OUTPUT_DIR/diamond_commands.txt"' not in script
        assert '"$OUTPUT_DIR/working_dir_path.txt"' not in script

    def test_passes_search_program_to_prepare(self, sample_config):
        # Prepare must invoke OrthoFinder with -S so the WorkingDirectory
        # is configured for the right search program. Resume relies on this.
        script = generate_prepare_script(sample_config)
        assert 'SEARCH_PROGRAM="diamond"' in script
        assert '-op -S "$SEARCH_PROGRAM"' in script

    def test_uses_strict_command_extraction_regex(self, sample_config):
        # The previous loose grep matched any line starting with "diamond"
        # — including the prose header "diamond commands that must be run"
        # — which broke the search array. The new regex must require a
        # subcommand token after "diamond".
        script = generate_prepare_script(sample_config)
        assert 'grep -E "^(diamond|blastp|makeblastdb)"' not in script
        assert "diamond[[:space:]]+(blastp|makedb)" in script
        assert "makeblastdb[[:space:]]" in script

    def test_validates_command_file_after_extraction(self, sample_config):
        # The prepare script must fail loudly if the command file is empty
        # or contains a non-command line, so the search array job never
        # receives garbage to `eval`.
        script = generate_prepare_script(sample_config)
        assert 'if [ ! -s "$COMMANDS_FILE" ]' in script
        assert "No valid DIAMOND/search commands were extracted" in script
        assert "BAD_LINES=" in script
        assert "Invalid lines found in command file" in script

    def test_strips_leading_whitespace_from_commands(self, sample_config):
        # OrthoFinder sometimes indents command lines. We strip leading
        # whitespace so the search loop can `eval` them directly.
        script = generate_prepare_script(sample_config)
        assert "sed 's/^[[:space:]]*//'" in script


class TestPrepareDbModeDetection:
    """Prepare must detect, on the user's cluster, whether this OrthoFinder
    build creates the search databases itself during -op or only emits the
    build commands, then persist that decision for later steps."""

    def test_detects_and_persists_db_mode(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert 'DB_MODE_FILE="$OUTPUT_PARENT/${RUN_NAME}_db_mode.txt"' in script
        assert 'echo "$OF_DB_MODE" > "$DB_MODE_FILE"' in script

    def test_uses_two_independent_signals(self, sample_config):
        # (1) build commands emitted by -op, and (2) DB files already present.
        script = generate_prepare_script(sample_config)
        assert "grep -E 'makedb|makeblastdb'" in script
        assert "grep -Ev 'makedb|makeblastdb'" in script
        assert "-name '*.dmnd'" in script

    def test_emit_build_commands_branch(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert 'OF_DB_MODE="emit_build_commands"' in script
        assert '[ "$DB_BUILD_COUNT" -gt 0 ]' in script

    def test_self_built_branch(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert 'OF_DB_MODE="self_built"' in script
        assert '[ "$DB_FILE_COUNT" -gt 0 ]' in script

    def test_fails_fast_when_mode_undeterminable(self, sample_config):
        # Neither build commands nor DB files found -> abort with a clear
        # diagnostic rather than guessing.
        script = generate_prepare_script(sample_config)
        assert (
            "Could not determine how OrthoFinder handled the search databases"
            in script
        )
        idx_err = script.index("Could not determine how OrthoFinder handled")
        assert "exit 1" in script[idx_err:]

    def test_detection_runs_after_workdir_located(self, sample_config):
        # The DB-file check needs $WORK_DIR, so detection must come after the
        # WorkingDirectory is located and its pointer written.
        script = generate_prepare_script(sample_config)
        assert script.index('echo "$WORK_DIR" > "$WORK_DIR_FILE"') < script.index(
            "OF_DB_MODE="
        )


class TestPrepareCompleteness:
    """Prepare verifies n databases (one per species) and n^2 searches, and
    builds the databases itself when OrthoFinder emitted the build commands."""

    def test_reads_species_count_from_speciesids(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert 'SPECIES_IDS_FILE="$WORK_DIR/SpeciesIDs.txt"' in script
        assert "grep -cE '^[0-9]+:'" in script

    def test_expected_searches_is_n_squared(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert "EXPECTED_SEARCHES=$(( N_SPECIES * N_SPECIES ))" in script

    def test_splits_into_search_and_db_build_files(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert "_search_commands.txt" in script
        assert "_db_build_commands.txt" in script
        assert "grep -Ev 'makedb|makeblastdb'" in script

    def test_fails_when_search_count_not_n_squared(self, sample_config):
        script = generate_prepare_script(sample_config)
        assert '[ "$SEARCH_COUNT" -ne "$EXPECTED_SEARCHES" ]' in script
        idx = script.index('[ "$SEARCH_COUNT" -ne "$EXPECTED_SEARCHES" ]')
        assert "exit 1" in script[idx:]

    def test_builds_dbs_in_emit_mode_and_verifies_count(self, sample_config):
        script = generate_prepare_script(sample_config)
        # exactly n build commands, each executed
        assert '[ "$DB_BUILD_COUNT" -ne "$N_SPECIES" ]' in script
        assert 'eval "$BUILD_CMD"' in script
        # after building, each species DB must be present (missing -> abort)
        assert '[ "$MISSING_DBS" -gt 0 ]' in script

    def test_verifies_each_species_db_exists_top_level(self, sample_config):
        # Per-species existence check at the top level; the dependencies/
        # self-test DB (a subdir) must be ignored. Regression for the first
        # cluster run that aborted on 115 vs 114 (find -maxdepth 2 swept up
        # WorkingDirectory/dependencies/diamondDBSpecies0.dmnd).
        script = generate_prepare_script(sample_config)
        assert "for SPECIES_ID in $(grep -oE '^[0-9]+'" in script
        assert "diamondDBSpecies${SPECIES_ID}.dmnd" in script
        assert '[ "$MISSING_DBS" -gt 0 ]' in script
        # detection counts DB files at maxdepth 1, never descending into
        # dependencies/.
        assert 'DB_FILE_COUNT=$(find "$WORK_DIR" -maxdepth 1' in script


class TestPrepareResourceDerivation:
    """Prepare CPU/memory/walltime are derived, not hardcoded: CPU is capped
    at min(16, cores); memory and walltime scale with proteome volume."""

    def test_cpus_capped_at_16(self):
        assert derive_prepare_cpus(44) == 16
        assert derive_prepare_cpus(16) == 16

    def test_cpus_uses_fewer_on_small_nodes(self):
        assert derive_prepare_cpus(12) == 12

    def test_memory_has_floor(self):
        assert derive_prepare_memory_mb(0, 0) == 4096

    def test_memory_scales_with_volume_and_rounds(self):
        small = derive_prepare_memory_mb(100 * 1024 * 1024, 20 * 1024 * 1024)
        large = derive_prepare_memory_mb(8000 * 1024 * 1024, 400 * 1024 * 1024)
        assert large > small
        assert large % 1024 == 0

    def test_walltime_is_clamped(self):
        assert derive_prepare_walltime(0) == "00:30:00"
        assert derive_prepare_walltime(10**12) == "04:00:00"

    def test_walltime_scales_with_volume(self):
        small = derive_prepare_walltime(200 * 1024 * 1024)
        large = derive_prepare_walltime(3000 * 1024 * 1024)
        assert large >= small

    def test_proteome_volume_reads_fasta_sizes(self, tmp_path):
        (tmp_path / "a.faa").write_text(">x\n" + "M" * 100 + "\n")
        (tmp_path / "b.fa").write_text(">y\n" + "M" * 300 + "\n")
        (tmp_path / "notes.txt").write_text("ignore me")
        total, largest, n = _proteome_volume(str(tmp_path))
        assert n == 2
        assert largest == (tmp_path / "b.fa").stat().st_size
        assert total == (
            (tmp_path / "a.faa").stat().st_size + (tmp_path / "b.fa").stat().st_size
        )

    def test_proteome_volume_missing_dir_returns_zero(self):
        assert _proteome_volume("/no/such/dir/xyz") == (0, 0, 0)

    def test_generator_embeds_derived_mem_and_time(self, tmp_path, sample_runtime):
        for i in range(4):
            (tmp_path / f"sp{i}.faa").write_text(">g\n" + "M" * 5000 + "\n")
        total, largest, _ = _proteome_volume(str(tmp_path))
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu", cpus_per_task=48),
            orthofinder=OrthoFinderConfig(
                input_dir=str(tmp_path),
                output_dir="/p/out/run",
                search_threads=16,
                analysis_threads=8,
            ),
            runtime=sample_runtime,
        )
        script = generate_prepare_script(cfg)
        assert f"--mem={derive_prepare_memory_mb(total, largest)}M" in script
        assert f"--time={derive_prepare_walltime(total)}" in script
        assert "--cpus-per-task=16" in script  # min(16, 48)

    def test_generator_falls_back_when_inputs_absent(self, sample_config):
        # sample_config.input_dir does not exist -> fallback mem/time.
        script = generate_prepare_script(sample_config)
        assert "--mem=4096M" in script
        assert "--time=02:00:00" in script


class TestCommandExtractionRegex:
    """End-to-end test that actually runs the prepare-script regex against
    synthetic OrthoFinder output. String-matching the script source can't
    catch regex bugs; running `grep` with the real pattern can."""

    SAMPLE_LOG = textwrap.dedent(
        """\
        OrthoFinder version 2.5.5
        ...analysing files...

        diamond commands that must be run
        diamond blastp --ignore-warnings -d db -q query.fa -o out.txt --more-sensitive
        diamond blastp --ignore-warnings -d db2 -q query2.fa -o out2.txt --more-sensitive
        diamond makedb --in species.fa -d diamondDBSpecies1
        makeblastdb -in species.fa -dbtype prot
            diamond blastp --ignore-warnings -d db3 -q query3.fa -o out3.txt
        some random text
        blastp commands that must be run
        Done!
        """
    )

    @pytest.mark.skipif(
        shutil.which("bash") is None or shutil.which("grep") is None,
        reason="bash and grep are required for extraction test",
    )
    def test_extraction_includes_real_commands_excludes_header(self, tmp_path):
        log = tmp_path / "prepare.log"
        log.write_text(self.SAMPLE_LOG)
        out = tmp_path / "commands.txt"

        # Use the exact extraction pipeline the prepare script uses.
        result = subprocess.run(
            [
                "bash",
                "-c",
                f"grep -E '{COMMAND_LINE_REGEX_EXTRACT}' \"$1\" "
                "| sed 's/^[[:space:]]*//' > \"$2\" || true",
                "_",
                str(log),
                str(out),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr

        lines = out.read_text().splitlines()

        # The OrthoFinder header lines must NOT appear. These are the
        # regressions — both the diamond and blastp variants of the
        # "X commands that must be run" prose.
        assert "diamond commands that must be run" not in lines
        assert "blastp commands that must be run" not in lines

        # Real commands must appear, with leading whitespace stripped.
        assert (
            "diamond blastp --ignore-warnings -d db -q query.fa -o out.txt "
            "--more-sensitive"
        ) in lines
        assert (
            "diamond blastp --ignore-warnings -d db2 -q query2.fa -o out2.txt "
            "--more-sensitive"
        ) in lines
        assert "diamond makedb --in species.fa -d diamondDBSpecies1" in lines
        assert "makeblastdb -in species.fa -dbtype prot" in lines
        # Indented line is captured and de-indented.
        assert (
            "diamond blastp --ignore-warnings -d db3 -q query3.fa -o out3.txt"
            in lines
        )
        # Prose and unrelated text are dropped.
        assert "some random text" not in lines

    @pytest.mark.parametrize(
        "bad_line",
        [
            "diamond commands that must be run",
            "blastp commands that must be run",
            "some random prose",
            "",
        ],
    )
    @pytest.mark.skipif(
        shutil.which("bash") is None,
        reason="bash is required to exercise the generated prepare script",
    )
    def test_prepare_script_validation_rejects_bad_line(
        self, sample_config, tmp_path, bad_line
    ):
        # Build a $COMMANDS_FILE containing the bad line that broke the
        # cluster (and a few siblings), then run the BAD_LINES validation
        # block the prepare script uses. It must exit non-zero.
        bad_commands = tmp_path / "diamond_commands.txt"
        bad_commands.write_text(bad_line + "\n")

        result = subprocess.run(
            [
                "bash",
                "-c",
                f"""
                set -uo pipefail
                COMMANDS_FILE="$1"
                BAD_LINES=$(grep -nEv '{COMMAND_LINE_REGEX}' "$COMMANDS_FILE" || true)
                if [ -n "$BAD_LINES" ]; then
                    echo "REJECTED"
                    exit 2
                fi
                echo "ACCEPTED"
                """,
                "_",
                str(bad_commands),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 2, (
            f"Expected validation to reject {bad_line!r}. "
            f"stdout={result.stdout!r} stderr={result.stderr!r}"
        )
        assert "REJECTED" in result.stdout

    @pytest.mark.skipif(
        shutil.which("bash") is None,
        reason="bash is required to exercise the generated prepare script",
    )
    def test_prepare_script_validation_accepts_real_commands(
        self, sample_config, tmp_path
    ):
        good_commands = tmp_path / "diamond_commands.txt"
        good_commands.write_text(
            "diamond blastp --ignore-warnings -d db -q q.fa -o o.txt\n"
            "diamond makedb --in s.fa -d dbS1\n"
            "makeblastdb -in s.fa -dbtype prot\n"
        )

        result = subprocess.run(
            [
                "bash",
                "-c",
                f"""
                set -uo pipefail
                COMMANDS_FILE="$1"
                BAD_LINES=$(grep -nEv '{COMMAND_LINE_REGEX}' "$COMMANDS_FILE" || true)
                if [ -n "$BAD_LINES" ]; then
                    echo "REJECTED: $BAD_LINES"
                    exit 2
                fi
                echo "ACCEPTED"
                """,
                "_",
                str(good_commands),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"Expected validation to accept real commands. "
            f"stdout={result.stdout!r} stderr={result.stderr!r}"
        )
        assert "ACCEPTED" in result.stdout


class TestSearchArrayStub:
    """The v1 search array has been removed; the v2 generator is a hard stub
    that raises NotImplementedError until the LPT/C-W-K-T sizing is built."""

    def test_raises_not_implemented(self, sample_config):
        with pytest.raises(NotImplementedError, match="search array"):
            generate_search_array_script(sample_config)


class TestResumeScript:
    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_resume" in generate_resume_script(sample_config)

    def test_uses_b_flag(self, sample_config):
        assert "-b " in generate_resume_script(sample_config)

    def test_uses_threads_from_config(self, sample_config):
        # sample_config has search_threads=16, analysis_threads=8
        script = generate_resume_script(sample_config)
        assert "TOTAL_THREADS=16" in script
        assert "ANALYSIS_THREADS=8" in script

    def test_does_not_pass_search_program(self, sample_config):
        # The resume command line should not contain "-S diamond" — search
        # program is already configured in WorkingDirectory by prepare phase.
        assert "-S diamond" not in generate_resume_script(sample_config)

    def test_finds_working_directory(self, sample_config):
        assert "WorkingDirectory" in generate_resume_script(sample_config)

    def test_reads_workdir_from_parent_not_output_dir(self, sample_config):
        # Must match the path prepare writes to.
        script = generate_resume_script(sample_config)
        assert (
            'WORK_DIR_FILE="$OUTPUT_PARENT/${RUN_NAME}_working_dir_path.txt"'
            in script
        )
        assert 'WORK_DIR_FILE="$OUTPUT_DIR/working_dir_path.txt"' not in script
        assert 'WORK_DIR="$(cat "$WORK_DIR_FILE")"' in script


class TestCrossCutting:
    def _all_scripts(self, config):
        # The search array is a v2 stub (raises), so cross-cutting checks cover
        # the two generators that actually produce scripts.
        return [
            generate_prepare_script(config),
            generate_resume_script(config),
        ]

    def test_all_scripts_have_set_euo_pipefail(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert "set -euo pipefail" in script

    def test_all_scripts_activate_conda(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert 'conda activate "$CONDA_ENV"' in script

    def test_all_scripts_no_conda_info_base_as_primary(self, sample_config):
        for script in self._all_scripts(sample_config):
            primary_end = script.index("elif command -v conda")
            assert "conda info --base" not in script[:primary_end]

    def test_all_scripts_have_export_all(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert "#SBATCH --export=ALL" in script

    def test_all_scripts_have_bootstrap_markers(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert "# ---- convgeno runtime bootstrap ----" in script
            assert "# ---- end convgeno runtime bootstrap ----" in script

    @pytest.mark.skipif(
        shutil.which("bash") is None,
        reason="bash not installed; cannot syntax-check generated scripts",
    )
    def test_all_scripts_are_syntactically_valid_bash(
        self, sample_config, tmp_path
    ):
        # Regression guard: every previous prepare-script bug was a bash
        # design error (mkdir/redirect order, wrong flags) that string-
        # matching tests missed. At minimum, every generated script must
        # parse under `bash -n`.
        names = ["prepare.sh", "resume.sh"]
        for name, script in zip(names, self._all_scripts(sample_config)):
            path = tmp_path / name
            path.write_text(script)
            result = subprocess.run(
                ["bash", "-n", str(path)],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, (
                f"Generated {name} has bash syntax errors:\n"
                f"{result.stderr}\n---script---\n{script}"
            )

    def test_raises_without_orthofinder_config(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu"),
            runtime=sample_runtime,
        )
        # prepare/resume validate the orthofinder block; the search array is a
        # v2 stub that raises NotImplementedError regardless of config.
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_prepare_script(config)
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_resume_script(config)
        with pytest.raises(NotImplementedError):
            generate_search_array_script(config)


class TestAccountOmissionInMultinodeScripts:
    @pytest.fixture()
    def no_account_config(self, sample_runtime) -> PipelineConfig:
        return PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                time_limit="72:00:00",
                account=None,
            ),
            orthofinder=OrthoFinderConfig(
                input_dir="/share/ceph/project/proteomes",
                output_dir="/share/ceph/project/results",
                search_threads=16,
                analysis_threads=8,
            ),
            runtime=sample_runtime,
        )

    def test_multinode_scripts_omit_null_account(self, no_account_config):
        prepare = generate_prepare_script(no_account_config)
        resume = generate_resume_script(no_account_config)
        for script in (prepare, resume):
            assert "--account" not in script

    def test_multinode_scripts_include_valid_account(self, sample_config):
        # sample_config has account="wym219"
        prepare = generate_prepare_script(sample_config)
        resume = generate_resume_script(sample_config)
        for script in (prepare, resume):
            assert "#SBATCH --account=wym219" in script

    @pytest.mark.parametrize("account", ["", "null", "None", "NULL"])
    def test_multinode_scripts_omit_stringy_null_accounts(
        self, account, no_account_config, sample_runtime
    ):
        config = PipelineConfig(
            project_dir=no_account_config.project_dir,
            conda_env=no_account_config.conda_env,
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                time_limit="72:00:00",
                account=account,
            ),
            orthofinder=no_account_config.orthofinder,
            runtime=sample_runtime,
        )
        prepare = generate_prepare_script(config)
        resume = generate_resume_script(config)
        for script in (prepare, resume):
            assert "--account" not in script


class TestResumeScriptTmpdirAndThreads:
    @pytest.fixture()
    def config_threads_unset(self, sample_runtime) -> PipelineConfig:
        return PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                time_limit="72:00:00",
            ),
            orthofinder=OrthoFinderConfig(
                input_dir="/share/ceph/project/proteomes",
                output_dir="/share/ceph/project/results",
                search_threads=16,
                analysis_threads=None,
            ),
            runtime=sample_runtime,
        )

    def test_resume_has_no_open_file_limit_logic(self, config_threads_unset):
        # The old raise-to-8192 / if-else ulimit block was removed; the fd
        # feasibility gate is reintroduced as a separate step.
        script = generate_resume_script(config_threads_unset)
        assert "ulimit -n" not in script
        assert "REQUESTED_OPEN_FILE_LIMIT" not in script
        assert "Open-file limit" not in script

    def test_resume_sets_tmpdir(self, config_threads_unset):
        script = generate_resume_script(config_threads_unset)
        assert 'export TMPDIR="$OUTPUT_PARENT/${RUN_NAME}_tmp"' in script
        assert 'mkdir -p "$TMPDIR"' in script

    def test_resume_cleans_tmpdir_on_success(self, config_threads_unset):
        script = generate_resume_script(config_threads_unset)
        assert 'rm -rf "$TMPDIR"' in script
        assert "Cleaning up TMPDIR" in script

    def test_resume_preserves_tmpdir_on_failure(self, config_threads_unset):
        script = generate_resume_script(config_threads_unset)
        assert "Preserving TMPDIR for debugging" in script

    def test_analysis_threads_default_to_1_when_unset(self, config_threads_unset):
        # analysis_threads=None → -a falls back to 1 (no auto-heuristic)
        script = generate_resume_script(config_threads_unset)
        assert "ANALYSIS_THREADS=1" in script

    def test_explicit_analysis_threads_override(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
            ),
            orthofinder=OrthoFinderConfig(
                input_dir="/in",
                output_dir="/out",
                search_threads=16,
                analysis_threads=6,
            ),
            runtime=sample_runtime,
        )
        script = generate_resume_script(config)
        assert "ANALYSIS_THREADS=6" in script

    def test_resume_uses_variable_not_hardcoded(self, config_threads_unset):
        script = generate_resume_script(config_threads_unset)
        assert '-a "$ANALYSIS_THREADS"' in script
        assert '-t "$TOTAL_THREADS"' in script

    def test_resume_exits_with_of_exit(self, config_threads_unset):
        script = generate_resume_script(config_threads_unset)
        assert 'exit "$OF_EXIT"' in script
