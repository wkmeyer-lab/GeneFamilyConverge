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
from convgeno.slurm.config import MultinodeConfig, PipelineConfig, SlurmConfig
from convgeno.slurm.multinode_generator import (
    COMMAND_LINE_REGEX,
    COMMAND_LINE_REGEX_EXTRACT,
    LPTResult,
    SearchArraySizing,
    SearchCommand,
    _parse_species_pair,
    _proteome_volume,
    build_all_pairs_search_commands,
    build_search_commands,
    compute_required_open_files,
    compute_search_array_sizing,
    derive_prepare_cpus,
    derive_prepare_memory_mb,
    derive_prepare_walltime,
    derive_search_array_sizing,
    derive_search_concurrency,
    derive_search_cpus,
    derive_search_cpus_from_discovery,
    derive_search_memory_mb,
    derive_search_task_count,
    derive_search_walltime,
    derive_search_waves,
    generate_prepare_script,
    generate_resume_script,
    generate_search_array_script,
    lpt_partition,
    read_species_fasta_sizes,
    read_species_ids,
    resolve_search_concurrency,
    resolve_search_task_count,
    search_task_manifest_name,
    within_task_concurrency,
    write_task_manifests,
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


class TestSearchArrayScript:
    """v2 search array: LPT-manifest-driven, C/p concurrency, idempotent skip,
    scratch (or direct-to-shared) outputs."""

    @pytest.fixture()
    def scratch_config(self, sample_runtime) -> PipelineConfig:
        return PipelineConfig(
            project_dir="/share/ceph/project",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                scratch_dir="/share/ceph/scratch",
                is_ephemeral_scratch=False,
            ),
            orthofinder=OrthoFinderConfig(
                input_dir="/share/ceph/project/proteomes",
                output_dir="/share/ceph/project/results",
                search_threads=16,
                analysis_threads=8,
            ),
            runtime=sample_runtime,
        )

    def test_shebang(self, sample_config):
        assert generate_search_array_script(sample_config).startswith("#!/bin/bash\n")

    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_search" in generate_search_array_script(
            sample_config
        )

    def test_array_header_matches_sizing(self, sample_config):
        sizing = derive_search_array_sizing(sample_config)
        throttle = min(sizing.concurrency, sizing.tasks)
        script = generate_search_array_script(sample_config)
        assert f"--array=0-{sizing.tasks - 1}%{throttle}" in script
        assert f"--cpus-per-task={sizing.cpus}" in script

    def test_reads_per_task_manifest(self, sample_config):
        script = generate_search_array_script(sample_config)
        assert 'search_task_${SLURM_ARRAY_TASK_ID}.txt' in script
        assert "_search_manifests" in script

    def test_idempotent_skip_with_gzip_test(self, sample_config):
        script = generate_search_array_script(sample_config)
        assert "gzip -t" in script
        assert "Blast([0-9]+)_([0-9]+)" in script  # pair parsed from Blast token
        assert "$STATUS_DIR/skipped" in script

    def test_concurrent_via_xargs(self, sample_config):
        script = generate_search_array_script(sample_config)
        assert "xargs" in script and "-P" in script
        assert "run_one" in script

    def test_uses_uo_pipefail_not_errexit(self, sample_config):
        # Per-command tolerance: the array must NOT abort on a failed command.
        script = generate_search_array_script(sample_config)
        assert "set -uo pipefail" in script
        assert "set -euo pipefail" not in script

    def test_direct_to_shared_when_no_scratch(self, sample_config):
        # sample_config has scratch_dir=None.
        script = generate_search_array_script(sample_config)
        assert 'OUTPUT_BASE="$WORK_DIR"' in script
        assert "Scratch not configured" in script
        assert "resolve_scratch_base" not in script

    def test_scratch_setup_when_configured(self, scratch_config):
        script = generate_search_array_script(scratch_config)
        assert "resolve_scratch_base" in script
        assert 'PREFERRED_SCRATCH="/share/ceph/scratch"' in script
        assert 'FALLBACK_SCRATCH="/tmp/scratch"' in script
        assert "salvage_scratch" in script
        assert 'rsync -a "$OUTPUT_BASE"/ "$WORK_DIR"/' in script

    def test_persistent_scratch_is_left(self, scratch_config):
        script = generate_search_array_script(scratch_config)
        assert "Persistent scratch left" in script

    def test_ephemeral_scratch_is_cleaned(self, sample_runtime):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(
                partition="hawkcpu",
                scratch_dir="/local/scratch",
                is_ephemeral_scratch=True,
            ),
            orthofinder=OrthoFinderConfig(input_dir="/in", output_dir="/out/run"),
            runtime=sample_runtime,
        )
        script = generate_search_array_script(cfg)
        assert "Cleaning up ephemeral scratch" in script

    def test_exits_nonzero_on_failures(self, sample_config):
        script = generate_search_array_script(sample_config)
        assert 'if [ "$FAILED" -gt 0 ]; then' in script


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

    def test_resume_passes_msa_tree_method(self, sample_config):
        # sample_config uses OrthoFinder defaults msa=mafft, tree=fasttree, so
        # the resume must select the MSA gene-tree method explicitly.
        script = generate_resume_script(sample_config)
        assert "-M msa -A mafft -T fasttree" in script

    def test_resume_method_matches_single_node(self, sample_config):
        # Scientific equivalence: the resume selects the SAME method/aligner/tree
        # as the single-node command. Single-node passes the aligner via a shell
        # variable (-A "$ORTHOFINDER_MSA_PROGRAM"), so compare effective values.
        from convgeno.slurm.script_generator import generate_orthofinder_script

        single = generate_orthofinder_script(sample_config)
        resume = generate_resume_script(sample_config)
        assert "-M msa -A mafft -T fasttree" in resume
        assert "-M msa" in single
        assert 'ORTHOFINDER_MSA_PROGRAM="mafft"' in single  # -A value
        assert "-T fasttree" in single

    def test_resume_omits_method_when_no_msa(self, sample_runtime):
        # msa_program unset -> OrthoFinder's dendroblast default (no -M/-A/-T),
        # matching the single-node behaviour for the same config.
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu", cpus_per_task=16),
            orthofinder=OrthoFinderConfig(
                input_dir="/in",
                output_dir="/out/run",
                search_threads=16,
                analysis_threads=8,
                msa_program="",
                tree_program="",
            ),
            runtime=sample_runtime,
        )
        script = generate_resume_script(cfg)
        assert "-M msa" not in script
        assert "-A " not in script
        assert "-T " not in script
        assert "dendroblast" in script  # the informational echo

    def test_resume_runs_completeness_gate_before_b(self, sample_config):
        # The gate must run, and must precede orthofinder -b so an incomplete
        # search set never reaches the analysis.
        script = generate_resume_script(sample_config)
        assert "verify_search_complete" in script
        gate = script.index("verify_search_complete")
        resume = script.index('orthofinder -b "$WORK_DIR"')
        assert gate < resume

    def test_resume_gate_aborts_without_b(self, sample_config):
        script = generate_resume_script(sample_config)
        assert "not running orthofinder -b" in script

    def test_resume_raises_ulimit_before_b(self, sample_config):
        # The fd gate must run (and set ulimit) before -b so orthofinder
        # inherits the raised soft limit.
        script = generate_resume_script(sample_config)
        assert 'ulimit -n "$REQUIRED_R"' in script
        assert "ulimit -Hn" in script
        fd = script.index("Open-file gate")
        resume = script.index('orthofinder -b "$WORK_DIR"')
        assert fd < resume

    def test_resume_fd_formula_matches_python(self, sample_config):
        # The bash required_r arithmetic must equal compute_required_open_files.
        script = generate_resume_script(sample_config)
        assert "(N_SPECIES * N_SPECIES * 11 + 9) / 10 + 1024" in script

    def test_resume_does_not_reduce_analysis_threads_for_fd(self, sample_config):
        # -a stays full (from config); the fd fix is ulimit, not fewer threads.
        script = generate_resume_script(sample_config)
        assert "ANALYSIS_THREADS=8" in script  # sample_config analysis_threads=8


class TestComputeRequiredOpenFiles:
    """required_r = ceil(n^2 * 1.1) + 1024, exact integer arithmetic."""

    def test_small_values(self):
        assert compute_required_open_files(0) == 1024  # 0 + headroom
        assert compute_required_open_files(1) == 1026  # ceil(1*1.1)=2, +1024

    def test_formula_matches_ceil(self):
        import math
        from fractions import Fraction

        # Exact rational reference — float ceil(n^2 * 1.1) mis-rounds when the
        # product lands on an integer (e.g. 100*1.1 == 110.00000000000001).
        for n in (2, 10, 114, 300):
            expected = math.ceil(Fraction(n * n) * Fraction(11, 10)) + 1024
            assert compute_required_open_files(n) == expected

    def test_covers_issue_571_figure(self):
        # Issue #571: 454 species need r ~= 206216. Our value must exceed it.
        assert compute_required_open_files(454) >= 206216

    def test_grows_quadratically(self):
        assert compute_required_open_files(200) > compute_required_open_files(100)

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="n must be"):
            compute_required_open_files(-1)

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
        return [
            generate_prepare_script(config),
            generate_search_array_script(config),
            generate_resume_script(config),
        ]

    def test_all_scripts_have_pipefail(self, sample_config):
        # prepare/resume abort on any error (set -euo pipefail); the search
        # array tolerates per-command failures (set -uo pipefail). Both share
        # the pipefail invariant.
        for script in self._all_scripts(sample_config):
            assert "pipefail" in script

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
        names = ["prepare.sh", "search.sh", "resume.sh"]
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
        # All three generators validate the orthofinder block first.
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_prepare_script(config)
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_resume_script(config)
        with pytest.raises(ValueError, match="OrthoFinder"):
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

    def test_resume_has_no_v1_open_file_block(self, config_threads_unset):
        # The v1 broken block (hardcoded raise-to-8192, -a tradeoff) must stay
        # gone; the v2 fd gate (ulimit -n required_r ~= n^2) replaces it.
        script = generate_resume_script(config_threads_unset)
        assert "8192" not in script
        assert "REQUESTED_OPEN_FILE_LIMIT" not in script

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


class TestLptPartition:
    """LPT cost-balanced bucketing of the n^2 search commands.

    Deterministic buckets from costs, full coverage with no duplicates, no
    empty tasks when T <= n^2, provable balance, and the two sizing byproducts:
    the biggest bucket's cost and the global-safe memory determinant.
    """

    @staticmethod
    def _cmds(specs: list[tuple[int, int, int]]) -> list[SearchCommand]:
        # specs: (line_index, query_bytes, db_bytes)
        return [SearchCommand(li, q, d) for li, q, d in specs]

    def test_known_small_case_is_exact(self):
        # Costs: 10, 8, 6, 5, 2 (total 31). Item li=3 carries the largest DB
        # (db=5) but a small cost, so it lands in the *smaller-cost* bucket --
        # this is the case that distinguishes global-safe mem from
        # biggest-bucket mem.
        cmds = self._cmds(
            [(0, 10, 1), (1, 8, 1), (2, 6, 1), (3, 1, 5), (4, 2, 1)]
        )
        res = lpt_partition(cmds, 2)
        assert isinstance(res, LPTResult)
        assert res.buckets == [[0, 3], [1, 2, 4]]
        assert res.bucket_costs == [15, 16]
        assert res.max_bucket_cost == 16

    def test_mem_determinant_is_global_not_biggest_bucket(self):
        # Largest DB (5) sits in bucket 0 (cost 15), NOT the biggest-cost
        # bucket (cost 16, all DBs = 1). Global-safe => 5, not 1.
        cmds = self._cmds(
            [(0, 10, 1), (1, 8, 1), (2, 6, 1), (3, 1, 5), (4, 2, 1)]
        )
        res = lpt_partition(cmds, 2)
        assert res.mem_determinant_bytes == 5

    def test_covers_every_command_exactly_once(self):
        cmds = self._cmds([(i, i + 1, (i % 3) + 1) for i in range(50)])
        res = lpt_partition(cmds, 7)
        flat = sorted(idx for bucket in res.buckets for idx in bucket)
        assert flat == list(range(50))

    def test_no_empty_bucket_when_buckets_le_commands(self):
        cmds = self._cmds([(i, i + 1, 1) for i in range(20)])
        res = lpt_partition(cmds, 20)
        assert all(len(bucket) >= 1 for bucket in res.buckets)

    def test_buckets_equal_commands_gives_one_each(self):
        cmds = self._cmds([(i, i + 1, 1) for i in range(6)])
        res = lpt_partition(cmds, 6)
        assert sorted(len(b) for b in res.buckets) == [1, 1, 1, 1, 1, 1]

    def test_more_buckets_than_commands_leaves_empties(self):
        cmds = self._cmds([(0, 3, 1), (1, 2, 1), (2, 1, 1)])
        res = lpt_partition(cmds, 5)
        assert len(res.buckets) == 5
        assert sum(1 for b in res.buckets if not b) == 2
        flat = sorted(idx for b in res.buckets for idx in b)
        assert flat == [0, 1, 2]

    def test_single_bucket_holds_everything(self):
        cmds = self._cmds([(0, 4, 1), (1, 3, 1), (2, 2, 1)])
        res = lpt_partition(cmds, 1)
        assert res.buckets == [[0, 1, 2]]
        assert res.bucket_costs == [9]
        assert res.max_bucket_cost == 9

    def test_bucket_costs_match_assigned_commands(self):
        cmds = self._cmds([(i, (i * 7) % 11 + 1, (i * 3) % 5 + 1) for i in range(40)])
        cost_by_index = {c.line_index: c.cost for c in cmds}
        res = lpt_partition(cmds, 6)
        for bucket, total in zip(res.buckets, res.bucket_costs):
            assert total == sum(cost_by_index[i] for i in bucket)
        assert res.max_bucket_cost == max(res.bucket_costs)

    def test_is_deterministic(self):
        cmds = self._cmds([(i, (i * 13) % 17 + 1, (i * 5) % 7 + 1) for i in range(60)])
        a = lpt_partition(cmds, 9)
        b = lpt_partition(list(reversed(cmds)), 9)
        assert a.buckets == b.buckets
        assert a.bucket_costs == b.bucket_costs

    def test_balance_within_greedy_bound(self):
        # Greedy least-loaded guarantees max_bucket <= total/T + max_item:
        # the last item added to the peak bucket found it at the minimum load,
        # which is <= the average <= total/T.
        cmds = self._cmds(
            [(i, (i * 31) % 97 + 1, (i * 17) % 53 + 1) for i in range(200)]
        )
        t = 11
        res = lpt_partition(cmds, t)
        total = sum(c.cost for c in cmds)
        max_item = max(c.cost for c in cmds)
        assert res.max_bucket_cost <= total / t + max_item

    def test_preserves_large_int_cost_precision(self):
        # 100 MB x 100 MB self-search = 1e16 > 2**53: must stay exact.
        big = 100 * 1024 * 1024
        cmds = self._cmds([(0, big, big), (1, 2, 1)])
        res = lpt_partition(cmds, 2)
        assert res.max_bucket_cost == big * big

    def test_raises_on_empty_commands(self):
        with pytest.raises(ValueError, match="non-empty"):
            lpt_partition([], 4)

    def test_raises_on_zero_buckets(self):
        cmds = self._cmds([(0, 1, 1)])
        with pytest.raises(ValueError, match="num_buckets"):
            lpt_partition(cmds, 0)


# A representative diamond blastp command as emitted into _search_commands.txt.
# The exact flags vary across OrthoFinder versions; the species pair is read
# from the version-stable Blast{i}_{j} output token.
def _diamond_cmd(i: int, j: int, workdir: str = "/w/WorkingDirectory") -> str:
    return (
        f"diamond blastp --ignore-warnings -d {workdir}/diamondDBSpecies{j} "
        f"-q {workdir}/Species{i}.fa -o {workdir}/Blast{i}_{j}.txt.gz "
        f"--more-sensitive -p 1 --quiet -e 0.001 --compress 1"
    )


class TestParseSpeciesPair:
    """Species-pair extraction from a single emitted search command."""

    def test_parses_query_i_and_db_j(self):
        # Blast{i}_{j}: i is the query species, j is the target database.
        assert _parse_species_pair(_diamond_cmd(7, 3)) == (7, 3)

    def test_db_species_substring_not_mistaken_for_query(self):
        # diamondDBSpecies3 contains "Species3" but must NOT be read as the
        # query; the query is Species7.fa.
        assert _parse_species_pair(_diamond_cmd(7, 3)) == (7, 3)

    def test_self_search_pair(self):
        assert _parse_species_pair(_diamond_cmd(4, 4)) == (4, 4)

    def test_blast_token_only_fallback(self):
        cmd = "diamond blastp -o /w/Blast2_5.txt.gz --opaque-future-flags"
        assert _parse_species_pair(cmd) == (2, 5)

    def test_query_db_tokens_when_no_blast_token(self):
        cmd = "blastp -query /w/Species8.fa -db /w/diamondDBSpecies1 -out foo"
        assert _parse_species_pair(cmd) == (8, 1)

    def test_query_disagrees_with_blast_raises(self):
        cmd = (
            "diamond blastp -d /w/diamondDBSpecies3 -q /w/Species8.fa "
            "-o /w/Blast7_3.txt.gz"
        )
        with pytest.raises(ValueError, match="disagrees"):
            _parse_species_pair(cmd)

    def test_db_disagrees_with_blast_raises(self):
        cmd = (
            "diamond blastp -d /w/diamondDBSpecies9 -q /w/Species7.fa "
            "-o /w/Blast7_3.txt.gz"
        )
        with pytest.raises(ValueError, match="disagrees"):
            _parse_species_pair(cmd)

    def test_unparseable_command_raises(self):
        with pytest.raises(ValueError, match="Could not determine"):
            _parse_species_pair("diamond blastp --help")


class TestBuildSearchCommands:
    """Assembling the SearchCommand list from command lines + proteome sizes."""

    def test_costs_use_query_and_db_sizes(self):
        # Asymmetric sizes prove query_bytes<-i, db_bytes<-j (the mem-critical
        # direction): Blast7_3 -> query=Species7 (200), db=Species3 (50).
        sizes = {3: 50, 7: 200}
        [cmd] = build_search_commands([_diamond_cmd(7, 3)], sizes)
        assert cmd.line_index == 0
        assert cmd.query_bytes == 200
        assert cmd.db_bytes == 50
        assert cmd.cost == 200 * 50

    def test_line_index_is_positional(self):
        sizes = {0: 10, 1: 20}
        lines = [_diamond_cmd(0, 1), _diamond_cmd(1, 0), _diamond_cmd(1, 1)]
        cmds = build_search_commands(lines, sizes)
        assert [c.line_index for c in cmds] == [0, 1, 2]

    def test_all_ordered_pairs_feed_lpt(self):
        # n=3 -> 9 ordered pairs; end-to-end into lpt_partition. The mem
        # determinant must equal the single largest proteome (200).
        sizes = {0: 100, 1: 200, 2: 50}
        lines = [_diamond_cmd(i, j) for i in range(3) for j in range(3)]
        cmds = build_search_commands(lines, sizes)
        assert len(cmds) == 9
        res = lpt_partition(cmds, 3)
        assert sum(len(b) for b in res.buckets) == 9
        assert res.mem_determinant_bytes == 200

    def test_blank_line_raises(self):
        with pytest.raises(ValueError, match="Blank line"):
            build_search_commands([_diamond_cmd(0, 0), "   "], {0: 10})

    def test_missing_species_size_raises(self):
        with pytest.raises(ValueError, match="No proteome size"):
            build_search_commands([_diamond_cmd(0, 5)], {0: 10})


class TestReadSpeciesFastaSizes:
    """Reading Species{id}.fa byte sizes from a WorkingDirectory."""

    def test_reads_sizes_and_ignores_other_files(self, tmp_path):
        (tmp_path / "Species0.fa").write_text(">a\n" + "M" * 100 + "\n")
        (tmp_path / "Species1.fa").write_text(">b\n" + "M" * 300 + "\n")
        # These must be ignored, not counted as proteomes.
        (tmp_path / "SpeciesIDs.txt").write_text("0: a.fa\n1: b.fa\n")
        (tmp_path / "diamondDBSpecies0.dmnd").write_bytes(b"\x00\x01\x02")
        sizes = read_species_fasta_sizes(tmp_path)
        assert set(sizes) == {0, 1}
        assert sizes[0] == (tmp_path / "Species0.fa").stat().st_size
        assert sizes[1] == (tmp_path / "Species1.fa").stat().st_size

    def test_no_species_files_raises(self, tmp_path):
        (tmp_path / "SpeciesIDs.txt").write_text("nothing useful")
        with pytest.raises(FileNotFoundError, match="No Species"):
            read_species_fasta_sizes(tmp_path)


class TestReadSpeciesIds:
    """Active species ids from SpeciesIDs.txt (commented species excluded)."""

    def test_reads_active_ids_sorted(self, tmp_path):
        (tmp_path / "SpeciesIDs.txt").write_text(
            "2: c.fa\n0: a.fa\n1: b.fa\n"
        )
        assert read_species_ids(tmp_path) == [0, 1, 2]

    def test_skips_commented_species(self, tmp_path):
        (tmp_path / "SpeciesIDs.txt").write_text(
            "0: a.fa\n#1: removed.fa\n2: c.fa\n"
        )
        assert read_species_ids(tmp_path) == [0, 2]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="SpeciesIDs.txt"):
            read_species_ids(tmp_path)

    def test_no_active_species_raises(self, tmp_path):
        (tmp_path / "SpeciesIDs.txt").write_text("#0: a.fa\n\n")
        with pytest.raises(ValueError, match="No active species"):
            read_species_ids(tmp_path)


class TestDeriveSearchCpus:
    """C = min(max_cpus // 3, min_cpus - 1), clamped >= 1, override wins."""

    def test_plan_example_heterogeneous(self):
        # biggest 52 / smallest 15 -> min(17, 14) = 14 (smallest-node cap).
        assert derive_search_cpus(52, 15) == 14

    def test_plan_example_uniform_16(self):
        # uniform 16-core -> min(5, 15) = 5 (third-of-biggest binds).
        assert derive_search_cpus(16, 16) == 5

    def test_smallest_node_cap_binds(self):
        # biggest 120 / smallest 8 -> min(40, 7) = 7.
        assert derive_search_cpus(120, 8) == 7

    def test_third_of_biggest_binds(self):
        # biggest 24 / smallest 100 -> min(8, 99) = 8.
        assert derive_search_cpus(24, 100) == 8

    def test_clamped_to_at_least_one(self):
        # Tiny nodes: min(0, 0) -> clamp to 1.
        assert derive_search_cpus(2, 1) == 1

    def test_empty_discovery_returns_one(self):
        assert derive_search_cpus(0, 0) == 1

    def test_override_wins(self):
        assert derive_search_cpus(52, 15, override=8) == 8

    def test_override_below_one_raises(self):
        with pytest.raises(ValueError, match="search_cpus override"):
            derive_search_cpus(52, 15, override=0)

    def test_from_discovery_dict(self):
        node = {"min_cpus_per_node": 15, "max_cpus_per_node": 52}
        assert derive_search_cpus_from_discovery(node) == 14

    def test_from_discovery_empty_dict(self):
        assert derive_search_cpus_from_discovery({}) == 1

    def test_from_discovery_override_passes_through(self):
        node = {"min_cpus_per_node": 15, "max_cpus_per_node": 52}
        assert derive_search_cpus_from_discovery(node, override=4) == 4


class TestDeriveSearchConcurrency:
    """W = the array %W throttle: default 8, override wins, QOS MaxJobs caps."""

    def test_default_is_eight(self):
        assert derive_search_concurrency() == 8

    def test_override_wins(self):
        assert derive_search_concurrency(override=12) == 12

    def test_override_below_one_raises(self):
        with pytest.raises(ValueError, match="array_throttle override"):
            derive_search_concurrency(override=0)

    def test_qos_caps_below_default(self):
        assert derive_search_concurrency(qos_max_jobs=3) == 3

    def test_qos_above_default_has_no_effect(self):
        assert derive_search_concurrency(qos_max_jobs=100) == 8

    def test_qos_caps_the_override(self):
        assert derive_search_concurrency(override=20, qos_max_jobs=5) == 5

    def test_qos_of_one_forces_serial(self):
        assert derive_search_concurrency(override=20, qos_max_jobs=1) == 1

    def test_qos_zero_means_unlimited(self):
        assert derive_search_concurrency(qos_max_jobs=0) == 8

    def test_qos_negative_ignored(self):
        assert derive_search_concurrency(qos_max_jobs=-5) == 8


class TestResolveSearchConcurrency:
    """W wired to discovery: derive_search_concurrency + detect_qos_max_jobs."""

    def test_discovered_qos_caps_default(self, monkeypatch):
        monkeypatch.setattr(
            "convgeno.slurm.multinode_generator.detect_qos_max_jobs",
            lambda partition: 3,
        )
        assert resolve_search_concurrency("hawkcpu") == 3

    def test_no_qos_falls_back_to_default(self, monkeypatch):
        monkeypatch.setattr(
            "convgeno.slurm.multinode_generator.detect_qos_max_jobs",
            lambda partition: None,
        )
        assert resolve_search_concurrency("hawkcpu") == 8

    def test_override_still_capped_by_qos(self, monkeypatch):
        monkeypatch.setattr(
            "convgeno.slurm.multinode_generator.detect_qos_max_jobs",
            lambda partition: 5,
        )
        assert resolve_search_concurrency("hawkcpu", override=20) == 5

    def test_override_used_when_qos_unknown(self, monkeypatch):
        monkeypatch.setattr(
            "convgeno.slurm.multinode_generator.detect_qos_max_jobs",
            lambda partition: None,
        )
        assert resolve_search_concurrency("hawkcpu", override=12) == 12


class TestDeriveSearchWaves:
    """K = waves: default 4, override wins."""

    def test_default_is_four(self):
        assert derive_search_waves() == 4

    def test_override_wins(self):
        assert derive_search_waves(override=6) == 6

    def test_override_below_one_raises(self):
        with pytest.raises(ValueError, match="waves .K. override"):
            derive_search_waves(override=0)


class TestDeriveSearchTaskCount:
    """T = min(K*W, num_commands, MaxArraySize), clamped >= 1."""

    def test_k_times_w_in_normal_regime(self):
        # 114 species -> n^2 = 12996; MaxArraySize 1001 -> T = K*W = 32.
        assert derive_search_task_count(12996, 8, 4, 1001) == 32

    def test_no_max_array_size_uses_num_commands_ceiling(self):
        assert derive_search_task_count(12996, 8, 4, None) == 32

    def test_capped_by_num_commands(self):
        # Tiny dataset: only 4 commands, so at most 4 tasks (< K*W = 32).
        assert derive_search_task_count(4, 8, 4, None) == 4

    def test_capped_by_max_array_size(self):
        assert derive_search_task_count(12996, 8, 4, 16) == 16

    def test_max_array_size_zero_ignored(self):
        assert derive_search_task_count(12996, 8, 4, 0) == 32

    def test_lower_bound_w_honored_when_feasible(self):
        # K=1 -> T = min(W, ceiling); with a large ceiling T == W (>= W floor).
        assert derive_search_task_count(12996, 8, 1, 1001) == 8

    def test_ceiling_below_w_wins_over_floor(self):
        # Only 4 commands but W=8: hard ceiling (4) beats the soft W floor.
        assert derive_search_task_count(4, 8, 4, None) == 4

    def test_clamped_to_at_least_one(self):
        assert derive_search_task_count(1, 1, 1, None) == 1

    def test_num_commands_below_one_raises(self):
        with pytest.raises(ValueError, match="num_commands"):
            derive_search_task_count(0, 8, 4)

    def test_concurrency_below_one_raises(self):
        with pytest.raises(ValueError, match="concurrency"):
            derive_search_task_count(100, 0, 4)

    def test_waves_below_one_raises(self):
        with pytest.raises(ValueError, match="waves"):
            derive_search_task_count(100, 8, 0)


class TestResolveSearchTaskCount:
    """T wired to discovery: derive_search_task_count + detect_max_array_size."""

    def test_uses_discovered_max_array_size(self, monkeypatch):
        monkeypatch.setattr(
            "convgeno.slurm.multinode_generator.detect_max_array_size",
            lambda: 16,
        )
        assert resolve_search_task_count(12996, 8, 4) == 16

    def test_no_limit_falls_back_to_num_commands(self, monkeypatch):
        monkeypatch.setattr(
            "convgeno.slurm.multinode_generator.detect_max_array_size",
            lambda: None,
        )
        assert resolve_search_task_count(12996, 8, 4) == 32


class TestWithinTaskConcurrency:
    """Concurrent diamond commands per task = C // p (>= 1)."""

    def test_p_one_equals_cpus(self):
        assert within_task_concurrency(14, 1) == 14

    def test_default_p_is_one(self):
        assert within_task_concurrency(14) == 14

    def test_integer_division(self):
        assert within_task_concurrency(16, 4) == 4
        assert within_task_concurrency(14, 3) == 4  # 14 // 3

    def test_floors_at_one(self):
        # p larger than C -> still one command at a time.
        assert within_task_concurrency(2, 8) == 1

    def test_invalid_cpus_raises(self):
        with pytest.raises(ValueError, match="cpus_per_task"):
            within_task_concurrency(0, 1)

    def test_invalid_threads_raises(self):
        with pytest.raises(ValueError, match="threads_per_command"):
            within_task_concurrency(8, 0)


class TestDeriveSearchWalltime:
    """bytes^2 bucket cost -> HH:MM:SS via C * throughput_const, clamped."""

    def test_none_when_not_estimable(self):
        assert derive_search_walltime(0, 14) is None
        assert derive_search_walltime(10**16, 0) is None
        assert derive_search_walltime(10**16, 14, throughput_const=0) is None

    def test_exact_mid_range_value(self):
        # 2.4e16 / (10 * 1e12) = 2400 s; * 1.5 margin = 3600 s = 01:00:00.
        got = derive_search_walltime(
            24_000_000_000_000_000, 10, throughput_const=1e12, margin=1.5
        )
        assert got == "01:00:00"

    def test_clamped_to_min(self):
        # Tiny bucket -> 15 min floor.
        assert derive_search_walltime(1, 14) == "00:15:00"

    def test_clamped_to_max(self):
        # Enormous bucket -> 72 h cap.
        assert derive_search_walltime(10**24, 1) == "72:00:00"

    def test_more_cpus_means_less_time(self):
        few = derive_search_walltime(10**17, 4, throughput_const=1e12)
        many = derive_search_walltime(10**17, 32, throughput_const=1e12)
        assert few > many  # HH:MM:SS strings compare correctly here

    def test_larger_margin_means_more_time(self):
        lo = derive_search_walltime(10**17, 8, throughput_const=1e12, margin=1.0)
        hi = derive_search_walltime(10**17, 8, throughput_const=1e12, margin=3.0)
        assert hi > lo


class TestDeriveSearchMemoryMb:
    """Largest target DB -> per-task MB; working-set floor, then DB scaling."""

    def test_none_when_not_estimable(self):
        assert derive_search_memory_mb(0, 14) is None
        assert derive_search_memory_mb(10 * 1024 * 1024, 0) is None

    def test_floor_dominates_for_small_proteomes(self):
        # 12 MB DB * 4 = 48 MB < 2048 floor -> per_command = 2048.
        # 14 * 2048 + 2048 base = 30720 MB.
        assert derive_search_memory_mb(12 * 1024 * 1024, 14) == 30720

    def test_db_scaling_for_large_databases(self):
        # 1 GB DB * 4 = 4096 MB > 2048 floor -> per_command = 4096.
        # 14 * 4096 + 2048 = 59392 MB.
        assert derive_search_memory_mb(1024 * 1024 * 1024, 14) == 59392

    def test_global_floor(self):
        # concurrency 1, tiny DB -> 1 * 2048 + 2048 = 4096 = floor.
        assert derive_search_memory_mb(1 * 1024 * 1024, 1) == 4096

    def test_rounds_up_to_gb(self):
        mem = derive_search_memory_mb(
            10 * 1024 * 1024, 7, base_mb=1500, per_command_floor_mb=1000
        )
        assert mem % 1024 == 0

    def test_scales_with_concurrency(self):
        small = derive_search_memory_mb(12 * 1024 * 1024, 4)
        large = derive_search_memory_mb(12 * 1024 * 1024, 32)
        assert large > small

    def test_overrides_apply(self):
        # db_safety binds: per_command = max(50*8, 100) = 400;
        # 12 * 400 + 1000 = 5800 -> round up to 6144 (above the 4 GB floor).
        mem = derive_search_memory_mb(
            50 * 1024 * 1024,
            12,
            db_safety=8.0,
            per_command_floor_mb=100,
            base_mb=1000,
        )
        assert mem == 6144


class TestWriteTaskManifests:
    """LPT buckets -> one self-contained command file per array task."""

    def test_name_matches_array_task_id(self):
        assert search_task_manifest_name(0) == "search_task_0.txt"
        assert search_task_manifest_name(31) == "search_task_31.txt"

    def test_writes_one_file_per_task_with_commands(self, tmp_path):
        cmds = [f"diamond blastp -o /w/Blast{i}_0.txt.gz -q x" for i in range(5)]
        buckets = [[0, 2, 4], [1, 3]]
        paths = write_task_manifests(buckets, cmds, tmp_path / "manifest")
        assert [p.name for p in paths] == ["search_task_0.txt", "search_task_1.txt"]
        assert paths[0].read_text().splitlines() == [cmds[0], cmds[2], cmds[4]]
        assert paths[1].read_text().splitlines() == [cmds[1], cmds[3]]

    def test_trailing_newline_for_clean_iteration(self, tmp_path):
        [path] = write_task_manifests([[0]], ["cmd-a"], tmp_path)
        assert path.read_text() == "cmd-a\n"

    def test_empty_bucket_writes_empty_file(self, tmp_path):
        paths = write_task_manifests([[0], []], ["only-cmd"], tmp_path)
        assert paths[1].exists()
        assert paths[1].read_text() == ""

    def test_creates_manifest_dir(self, tmp_path):
        target = tmp_path / "nested" / "manifest"
        write_task_manifests([[0]], ["cmd"], target)
        assert target.is_dir()

    def test_every_command_appears_exactly_once(self, tmp_path):
        cmds = [f"cmd{i}" for i in range(20)]
        # A plausible LPT-style partition of all 20 indices over 6 tasks.
        buckets = [
            [0, 6, 12, 18], [1, 7, 13, 19], [2, 8, 14],
            [3, 9, 15], [4, 10, 16], [5, 11, 17],
        ]
        paths = write_task_manifests(buckets, cmds, tmp_path)
        assert len(paths) == 6
        seen = [line for p in paths for line in p.read_text().splitlines()]
        assert sorted(seen) == sorted(cmds)

    def test_out_of_range_index_raises(self, tmp_path):
        with pytest.raises(ValueError, match="out of range"):
            write_task_manifests([[0, 5]], ["only-one-cmd"], tmp_path)

    def test_end_to_end_from_lpt(self, tmp_path):
        # build_search_commands -> lpt_partition -> write_task_manifests, and
        # confirm the manifests together reproduce exactly the input commands.
        sizes = {0: 100, 1: 200, 2: 50}
        lines = [
            f"diamond blastp -d /w/diamondDBSpecies{j} -q /w/Species{i}.fa "
            f"-o /w/Blast{i}_{j}.txt.gz -p 1"
            for i in range(3)
            for j in range(3)
        ]
        cmds = build_search_commands(lines, sizes)
        result = lpt_partition(cmds, 3)
        paths = write_task_manifests(result.buckets, lines, tmp_path)
        assert len(paths) == 3
        seen = [line for p in paths for line in p.read_text().splitlines()]
        assert sorted(seen) == sorted(lines)


class TestBuildAllPairsSearchCommands:
    """Generate-time n^2 estimate commands from a proteome-size list."""

    def test_count_and_pairs(self):
        cmds = build_all_pairs_search_commands([100, 200, 50])
        assert len(cmds) == 9  # n^2
        assert [c.line_index for c in cmds] == list(range(9))
        # cost is |S_i| * |S_j| over the full cartesian product.
        assert max(c.cost for c in cmds) == 200 * 200
        assert {c.db_bytes for c in cmds} == {100, 200, 50}

    def test_empty_sizes(self):
        assert build_all_pairs_search_commands([]) == []


class TestComputeSearchArraySizing:
    """The whole C/W/K/T + time/mem bundle from sizes + cluster numbers."""

    def test_normal_regime(self):
        sizing = compute_search_array_sizing(
            [10_000_000] * 114,  # 114 ~10 MB proteomes
            max_cpus_per_node=52,
            min_cpus_per_node=15,
            qos_max_jobs=None,
            max_array_size=1001,
            fallback_time="72:00:00",
            fallback_mem="16000M",
        )
        assert sizing.cpus == 14  # min(52//3, 15-1)
        assert sizing.concurrency == 8
        assert sizing.waves == 4
        assert sizing.tasks == 32  # K*W, well under n^2=12996 and MaxArraySize
        assert sizing.within_task_parallel == 14  # C/p, p=1
        assert sizing.mem.endswith("M")

    def test_qos_caps_concurrency(self):
        sizing = compute_search_array_sizing(
            [10_000_000] * 20,
            max_cpus_per_node=52,
            min_cpus_per_node=15,
            qos_max_jobs=3,
            max_array_size=1001,
            fallback_time="72:00:00",
            fallback_mem="16000M",
        )
        assert sizing.concurrency == 3  # W capped by QOS

    def test_off_cluster_falls_back(self):
        # No sizes, no discovery: safe fallbacks, still a valid array.
        sizing = compute_search_array_sizing(
            [],
            max_cpus_per_node=0,
            min_cpus_per_node=0,
            qos_max_jobs=None,
            max_array_size=None,
            fallback_time="72:00:00",
            fallback_mem="16000M",
        )
        assert sizing.cpus == 1
        assert sizing.concurrency == 8
        assert sizing.tasks == 32  # K*W (no n^2 cap available)
        assert sizing.time_limit == "72:00:00"  # fallback
        assert sizing.mem == "16000M"  # fallback

    def test_overrides_win(self):
        sizing = compute_search_array_sizing(
            [10_000_000] * 50,
            max_cpus_per_node=52,
            min_cpus_per_node=15,
            qos_max_jobs=None,
            max_array_size=1001,
            fallback_time="72:00:00",
            fallback_mem="16000M",
            overrides=MultinodeConfig(
                search_cpus=8,
                array_throttle=4,
                waves=2,
                search_time_limit="03:00:00",
                search_mem="24000M",
            ),
        )
        assert sizing.cpus == 8
        assert sizing.concurrency == 4
        assert sizing.waves == 2
        assert sizing.tasks == 8  # K*W = 2*4
        assert sizing.time_limit == "03:00:00"
        assert sizing.mem == "24000M"

    def test_returns_dataclass(self):
        sizing = compute_search_array_sizing(
            [], 0, 0, None, None, fallback_time="1:00:00", fallback_mem="8000M"
        )
        assert isinstance(sizing, SearchArraySizing)
