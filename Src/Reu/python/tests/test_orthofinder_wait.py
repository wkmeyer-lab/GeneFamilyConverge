"""Tests for the --wait job-polling logic in convgeno.cli.orthofinder_cmd."""

from __future__ import annotations

import subprocess
from unittest.mock import patch

from convgeno.cli.orthofinder_cmd import (
    _normalize_state,
    sacct_job_state,
    wait_for_completion,
)

NOOP_SLEEP = lambda _s: None  # noqa: E731


class TestNormalizeState:
    def test_plain(self):
        assert _normalize_state("COMPLETED") == "COMPLETED"

    def test_cancelled_by_uid(self):
        assert _normalize_state("CANCELLED by 12345") == "CANCELLED"

    def test_trailing_plus(self):
        assert _normalize_state("COMPLETED+") == "COMPLETED"

    def test_whitespace(self):
        assert _normalize_state("  RUNNING  ") == "RUNNING"


class TestWaitForCompletion:
    def test_completes(self):
        states = iter(["PENDING", "RUNNING", "COMPLETED"])
        rc = wait_for_completion(
            "1", sleep=NOOP_SLEEP, state_fn=lambda _j: next(states)
        )
        assert rc == 0

    def test_fails(self):
        states = iter(["RUNNING", "FAILED"])
        rc = wait_for_completion(
            "1", sleep=NOOP_SLEEP, state_fn=lambda _j: next(states)
        )
        assert rc == 1

    def test_dead_dependency_is_failure(self):
        # Upstream job failed -> resume is CANCELLED and never runs. Must not hang.
        rc = wait_for_completion("1", sleep=NOOP_SLEEP, state_fn=lambda _j: "CANCELLED")
        assert rc == 1

    def test_unknown_gives_up(self):
        rc = wait_for_completion(
            "1", max_unknown_polls=3, sleep=NOOP_SLEEP, state_fn=lambda _j: None
        )
        assert rc == 1


def _fake_completed(stdout: str = "", returncode: int = 0):
    return subprocess.CompletedProcess(
        args=[], returncode=returncode, stdout=stdout, stderr=""
    )


class TestSacctJobState:
    def test_completed(self):
        with patch(
            "convgeno.cli.orthofinder_cmd.subprocess.run",
            return_value=_fake_completed("COMPLETED\n"),
        ):
            assert sacct_job_state("1") == "COMPLETED"

    def test_array_active_wins(self):
        with patch(
            "convgeno.cli.orthofinder_cmd.subprocess.run",
            return_value=_fake_completed("COMPLETED\nRUNNING\n"),
        ):
            assert sacct_job_state("1") == "RUNNING"

    def test_array_failure_detected(self):
        with patch(
            "convgeno.cli.orthofinder_cmd.subprocess.run",
            return_value=_fake_completed("COMPLETED\nFAILED\n"),
        ):
            assert sacct_job_state("1") == "FAILED"

    def test_empty_is_unknown(self):
        with patch(
            "convgeno.cli.orthofinder_cmd.subprocess.run",
            return_value=_fake_completed("\n"),
        ):
            assert sacct_job_state("1") is None

    def test_nonzero_is_unknown(self):
        with patch(
            "convgeno.cli.orthofinder_cmd.subprocess.run",
            return_value=_fake_completed("", returncode=1),
        ):
            assert sacct_job_state("1") is None

    def test_missing_sacct_is_unknown(self):
        with patch(
            "convgeno.cli.orthofinder_cmd.subprocess.run",
            side_effect=FileNotFoundError(),
        ):
            assert sacct_job_state("1") is None
