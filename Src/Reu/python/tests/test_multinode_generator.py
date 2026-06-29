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
    choose_analysis_threads,
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
            time_limit="48:00:00",
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
        assert "--cpus-per-task=4" in script

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


def _search_script(config, **kwargs):
    """Helper: call generate_search_array_script with sensible defaults."""
    kwargs.setdefault("commands_per_task", 50)
    kwargs.setdefault("array_max", 99)
    return generate_search_array_script(config, **kwargs)


class TestSearchArrayScript:
    def test_is_array_job(self, sample_config):
        assert "--array=" in _search_script(sample_config)

    def test_job_name(self, sample_config):
        assert "--job-name=convgeno_of_search" in _search_script(sample_config)

    def test_reads_commands_file(self, sample_config):
        assert "diamond_commands.txt" in _search_script(sample_config)

    def test_reads_commands_from_parent_not_output_dir(self, sample_config):
        # Must match the path prepare writes to. The previous bug had
        # both prepare and search referencing $OUTPUT_DIR/diamond_commands.txt
        # which fails on the prepare side.
        script = _search_script(sample_config)
        assert (
            'COMMANDS_FILE="$OUTPUT_PARENT/${RUN_NAME}_diamond_commands.txt"'
            in script
        )
        assert 'COMMANDS_FILE="$OUTPUT_DIR/diamond_commands.txt"' not in script

    def test_uses_slurm_array_task_id(self, sample_config):
        assert "SLURM_ARRAY_TASK_ID" in _search_script(sample_config)

    def test_exits_if_no_work(self, sample_config):
        assert "exit 0" in _search_script(sample_config)

    def test_counts_failures(self, sample_config):
        assert "FAILED" in _search_script(sample_config)

    def test_uses_search_threads(self, sample_config):
        assert "--cpus-per-task=16" in _search_script(sample_config)

    def test_custom_commands_per_task(self, sample_config):
        script = _search_script(sample_config, commands_per_task=100)
        assert "COMMANDS_PER_TASK=100" in script

    def test_uses_computed_array_max(self, sample_config):
        # The hardcoded 0-9999 was rejected by SLURM clusters with
        # "Invalid job array specification". The array range must come
        # from caller-supplied array_max.
        script = generate_search_array_script(
            sample_config, commands_per_task=50, array_max=259
        )
        assert "#SBATCH --array=0-259" in script
        assert "#SBATCH --array=0-9999" not in script

    def test_array_max_zero_is_valid_single_task(self, sample_config):
        script = generate_search_array_script(
            sample_config, commands_per_task=50, array_max=0
        )
        assert "#SBATCH --array=0-0" in script

    def test_invalid_commands_per_task_raises(self, sample_config):
        with pytest.raises(ValueError, match="commands_per_task"):
            generate_search_array_script(
                sample_config, commands_per_task=0, array_max=259
            )

    def test_invalid_array_max_raises(self, sample_config):
        with pytest.raises(ValueError, match="array_max"):
            generate_search_array_script(
                sample_config, commands_per_task=50, array_max=-1
            )

    def test_missing_array_max_raises(self, sample_config):
        with pytest.raises(ValueError, match="array_max"):
            generate_search_array_script(
                sample_config, commands_per_task=50, array_max=None
            )


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
        return [
            generate_prepare_script(config),
            _search_script(config),
            generate_resume_script(config),
        ]

    def test_all_scripts_have_set_euo_pipefail(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert "set -euo pipefail" in script

    def test_all_scripts_activate_conda(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert 'conda activate "/home/prm526/.conda/envs/convgeno"' in script

    def test_all_scripts_no_conda_info_base(self, sample_config):
        for script in self._all_scripts(sample_config):
            assert "conda info --base" not in script

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
        # The orthofinder-missing check fires before array_max validation
        # in generate_search_array_script, so all three raise the same
        # "OrthoFinder" error even if we pass valid array_max.
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_prepare_script(config)
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_search_array_script(
                config, commands_per_task=50, array_max=99
            )
        with pytest.raises(ValueError, match="OrthoFinder"):
            generate_resume_script(config)


class TestAccountOmissionInMultinodeScripts:
    @pytest.fixture()
    def no_account_config(self, sample_runtime) -> PipelineConfig:
        return PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                time_limit="48:00:00",
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
        search = _search_script(no_account_config)
        resume = generate_resume_script(no_account_config)
        for script in (prepare, search, resume):
            assert "--account" not in script

    def test_multinode_scripts_include_valid_account(self, sample_config):
        # sample_config has account="wym219"
        prepare = generate_prepare_script(sample_config)
        search = _search_script(sample_config)
        resume = generate_resume_script(sample_config)
        for script in (prepare, search, resume):
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
                time_limit="48:00:00",
                account=account,
            ),
            orthofinder=no_account_config.orthofinder,
            runtime=sample_runtime,
        )
        prepare = generate_prepare_script(config)
        search = _search_script(config)
        resume = generate_resume_script(config)
        for script in (prepare, search, resume):
            assert "--account" not in script


class TestChooseAnalysisThreads:
    def test_none_limit_returns_conservative(self):
        assert choose_analysis_threads(16, None) == 2

    def test_low_limit_returns_1(self):
        assert choose_analysis_threads(16, 2048) == 1

    def test_medium_limit_returns_2(self):
        assert choose_analysis_threads(16, 4096) == 2

    def test_high_limit_returns_4(self):
        assert choose_analysis_threads(16, 8192) == 4

    def test_very_high_limit_returns_8(self):
        assert choose_analysis_threads(16, 16384) == 8

    def test_capped_by_cpus(self):
        assert choose_analysis_threads(2, 65536) == 2

    def test_single_cpu(self):
        assert choose_analysis_threads(1, 65536) == 1


class TestResumeScriptUlimitAndTmpdir:
    @pytest.fixture()
    def config_with_limit(self, sample_runtime) -> PipelineConfig:
        return PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                time_limit="48:00:00",
                open_file_limit=8192,
            ),
            orthofinder=OrthoFinderConfig(
                input_dir="/share/ceph/project/proteomes",
                output_dir="/share/ceph/project/results",
                search_threads=16,
                analysis_threads=None,
            ),
            runtime=sample_runtime,
        )

    @pytest.fixture()
    def config_no_limit(self, sample_runtime) -> PipelineConfig:
        return PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                time_limit="48:00:00",
                open_file_limit=None,
            ),
            orthofinder=OrthoFinderConfig(
                input_dir="/share/ceph/project/proteomes",
                output_dir="/share/ceph/project/results",
                search_threads=16,
                analysis_threads=None,
            ),
            runtime=sample_runtime,
        )

    def test_resume_contains_ulimit_n(self, config_with_limit):
        script = generate_resume_script(config_with_limit)
        assert "ulimit -n" in script

    def test_ulimit_wrapped_in_if_not_bare(self, config_with_limit):
        # The ulimit call must be safe under set -e. It should be
        # inside an if block so a refusal doesn't abort the script.
        script = generate_resume_script(config_with_limit)
        assert 'if ulimit -n "$REQUESTED_OPEN_FILE_LIMIT"' in script

    def test_ulimit_fallback_warning(self, config_with_limit):
        script = generate_resume_script(config_with_limit)
        assert "Could not raise open-file limit" in script

    def test_no_limit_prints_current(self, config_no_limit):
        script = generate_resume_script(config_no_limit)
        assert "Open-file limit (current):" in script
        assert "REQUESTED_OPEN_FILE_LIMIT" not in script

    def test_resume_sets_tmpdir(self, config_with_limit):
        script = generate_resume_script(config_with_limit)
        assert 'export TMPDIR="$OUTPUT_PARENT/${RUN_NAME}_tmp"' in script
        assert 'mkdir -p "$TMPDIR"' in script

    def test_resume_cleans_tmpdir_on_success(self, config_with_limit):
        script = generate_resume_script(config_with_limit)
        assert 'rm -rf "$TMPDIR"' in script
        assert "Cleaning up TMPDIR" in script

    def test_resume_preserves_tmpdir_on_failure(self, config_with_limit):
        script = generate_resume_script(config_with_limit)
        assert "Preserving TMPDIR for debugging" in script

    def test_auto_analysis_threads_with_limit_8192(self, config_with_limit):
        # open_file_limit=8192, analysis_threads=None → auto = 4
        script = generate_resume_script(config_with_limit)
        assert "ANALYSIS_THREADS=4" in script
        assert "-a 8" not in script

    def test_auto_analysis_threads_without_limit(self, config_no_limit):
        # open_file_limit=None, analysis_threads=None → auto = 2
        script = generate_resume_script(config_no_limit)
        assert "ANALYSIS_THREADS=2" in script

    def test_explicit_analysis_threads_override(self, sample_runtime):
        config = PipelineConfig(
            project_dir="/share/ceph/project",
            conda_env="convgeno",
            slurm=SlurmConfig(
                partition="hawkcpu",
                cpus_per_task=16,
                open_file_limit=8192,
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

    def test_resume_uses_variable_not_hardcoded(self, config_with_limit):
        script = generate_resume_script(config_with_limit)
        assert '-a "$ANALYSIS_THREADS"' in script
        assert '-t "$TOTAL_THREADS"' in script

    def test_resume_exits_with_of_exit(self, config_with_limit):
        script = generate_resume_script(config_with_limit)
        assert 'exit "$OF_EXIT"' in script
