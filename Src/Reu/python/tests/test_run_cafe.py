"""Tests for the Src/Loc/scripts/run_cafe.py CLI (CAFE-5 execution step).

The script is a standalone Loc entry point (not part of the ``convgeno`` package),
so it is loaded from its file path via importlib. CAFE-5 itself is not installed
in test environments, so the executable check and the command runner are mocked.
"""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path
from unittest.mock import patch

_SCRIPT = (
    Path(__file__).resolve().parents[4] / "Src" / "Loc" / "scripts" / "run_cafe.py"
)
_spec = importlib.util.spec_from_file_location("run_cafe", _SCRIPT)
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)


def _fake_run(write_results=True):
    """Return a fake command_runner.run that (optionally) writes CAFE results."""

    def _run(cmd, log_dir=None, **kwargs):
        if write_results:
            out = Path(cmd[cmd.index("-o") + 1])
            out.mkdir(parents=True, exist_ok=True)
            (out / "Base_results.txt").write_text("ok", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    return _run


def _inputs(tmp_path: Path):
    counts = tmp_path / "cafe_input.tsv"
    tree = tmp_path / "species_tree.nwk"
    counts.write_text("Desc\tFamily ID\tsp1\n", encoding="utf-8")
    tree.write_text("(sp1:1,sp2:1);\n", encoding="utf-8")
    return counts, tree, tmp_path / "cafe_results"


class TestRunCafe:
    def test_missing_cafe_binary_errors(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        with patch.object(rc, "check_tool_available", return_value=False):
            assert rc.run_cafe(counts, tree, out) == 1

    def test_runs_and_builds_command(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        with (
            patch.object(rc, "check_tool_available", return_value=True),
            patch.object(rc.cafe, "validate_input", return_value=[]),
            patch.object(rc, "run_command", side_effect=_fake_run()) as mock_run,
        ):
            code = rc.run_cafe(counts, tree, out)
        assert code == 0
        cmd = mock_run.call_args.args[0]
        assert "-i" in cmd and "-t" in cmd
        assert cmd[cmd.index("-o") + 1] == str(out)
        assert "-y" not in cmd  # no phenotype tree given

    def test_input_validation_failure_errors(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        with (
            patch.object(rc, "check_tool_available", return_value=True),
            patch.object(rc.cafe, "validate_input", return_value=["bad header"]),
            patch.object(rc, "run_command", side_effect=_fake_run()) as mock_run,
        ):
            assert rc.run_cafe(counts, tree, out) == 1
        mock_run.assert_not_called()

    def test_missing_lambda_tree_drops_y(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        missing_lambda = tmp_path / "lambda_tree.nwk"  # never created
        with (
            patch.object(rc, "check_tool_available", return_value=True),
            patch.object(rc.cafe, "validate_input", return_value=[]),
            patch.object(rc, "run_command", side_effect=_fake_run()) as mock_run,
        ):
            code = rc.run_cafe(counts, tree, out, lambda_tree=missing_lambda)
        assert code == 0
        assert "-y" not in mock_run.call_args.args[0]

    def test_present_but_invalid_lambda_tree_errors(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        lam = tmp_path / "lambda_tree.nwk"
        lam.write_text("(sp1:1,sp2:1);\n", encoding="utf-8")
        with (
            patch.object(rc, "check_tool_available", return_value=True),
            patch.object(rc.cafe, "validate_input", return_value=[]),
            patch.object(rc.cafe, "validate_lambda_tree", return_value=["mismatch"]),
            patch.object(rc, "run_command", side_effect=_fake_run()) as mock_run,
        ):
            assert rc.run_cafe(counts, tree, out, lambda_tree=lam) == 1
        mock_run.assert_not_called()

    def test_valid_lambda_tree_adds_y(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        lam = tmp_path / "lambda_tree.nwk"
        lam.write_text("(sp1:1,sp2:1);\n", encoding="utf-8")
        with (
            patch.object(rc, "check_tool_available", return_value=True),
            patch.object(rc.cafe, "validate_input", return_value=[]),
            patch.object(rc.cafe, "validate_lambda_tree", return_value=[]),
            patch.object(rc, "run_command", side_effect=_fake_run()) as mock_run,
        ):
            code = rc.run_cafe(counts, tree, out, lambda_tree=lam)
        assert code == 0
        cmd = mock_run.call_args.args[0]
        assert "-y" in cmd and cmd[cmd.index("-y") + 1] == str(lam)

    def test_output_validation_failure_errors(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        with (
            patch.object(rc, "check_tool_available", return_value=True),
            patch.object(rc.cafe, "validate_input", return_value=[]),
            # run "succeeds" but writes no *_results.txt -> output check fails.
            patch.object(rc, "run_command", side_effect=_fake_run(write_results=False)),
        ):
            assert rc.run_cafe(counts, tree, out) == 1

    def test_gamma_and_poisson_flags(self, tmp_path: Path):
        counts, tree, out = _inputs(tmp_path)
        with (
            patch.object(rc, "check_tool_available", return_value=True),
            patch.object(rc.cafe, "validate_input", return_value=[]),
            patch.object(rc, "run_command", side_effect=_fake_run()) as mock_run,
        ):
            code = rc.run_cafe(counts, tree, out, gamma_categories=3, poisson=True)
        assert code == 0
        cmd = mock_run.call_args.args[0]
        assert "-k" in cmd and cmd[cmd.index("-k") + 1] == "3"
        assert "-p" in cmd
