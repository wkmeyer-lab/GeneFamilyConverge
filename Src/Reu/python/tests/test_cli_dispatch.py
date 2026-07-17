"""Tests for the convgeno CLI dispatch + execution-mode resolution.

Multi-node is the default; single-node is opt-in via --mode single-node or the
--single-node alias.
"""

from __future__ import annotations

import argparse
import sys

import pytest

from convgeno.cli import _resolve_orthofinder_mode, main


def _args(mode=None, multinode=False, single_node=False) -> argparse.Namespace:
    return argparse.Namespace(mode=mode, multinode=multinode, single_node=single_node)


class TestResolveMode:
    def setup_method(self):
        self.parser = argparse.ArgumentParser()

    def test_default_is_multinode(self):
        assert _resolve_orthofinder_mode(self.parser, _args()) == "multinode"

    def test_multinode_alias(self):
        assert _resolve_orthofinder_mode(
            self.parser, _args(multinode=True)
        ) == "multinode"

    def test_mode_multinode_explicit(self):
        assert _resolve_orthofinder_mode(
            self.parser, _args(mode="multinode")
        ) == "multinode"

    def test_single_node_alias(self):
        assert _resolve_orthofinder_mode(
            self.parser, _args(single_node=True)
        ) == "single-node"

    def test_mode_single_node_explicit(self):
        assert _resolve_orthofinder_mode(
            self.parser, _args(mode="single-node")
        ) == "single-node"

    def test_conflict_multinode_and_single_node_errors(self):
        with pytest.raises(SystemExit):
            _resolve_orthofinder_mode(
                self.parser, _args(multinode=True, single_node=True)
            )

    def test_conflict_mode_single_with_multinode_flag_errors(self):
        with pytest.raises(SystemExit):
            _resolve_orthofinder_mode(
                self.parser, _args(mode="single-node", multinode=True)
            )


class TestMainDispatch:
    """End-to-end: argv -> mode resolution -> the correct handler is invoked."""

    def _run(self, monkeypatch, argv):
        monkeypatch.setattr(sys, "argv", ["convgeno", *argv])
        called = {}
        import convgeno.cli.orthofinder_cmd as oc

        for name in (
            "run_generate",
            "run_generate_multinode",
            "run_submit",
            "run_submit_multinode",
        ):
            monkeypatch.setattr(
                oc, name, lambda *a, _n=name, **k: called.setdefault(_n, (a, k))
            )
        main()
        return called

    def test_generate_defaults_to_multinode(self, monkeypatch):
        called = self._run(monkeypatch, ["orthofinder", "generate", "--config", "c"])
        assert set(called) == {"run_generate_multinode"}

    def test_generate_single_node_alias(self, monkeypatch):
        called = self._run(
            monkeypatch, ["orthofinder", "generate", "--single-node", "--config", "c"]
        )
        assert set(called) == {"run_generate"}

    def test_generate_mode_single_node(self, monkeypatch):
        called = self._run(
            monkeypatch,
            ["orthofinder", "generate", "--mode", "single-node", "--config", "c"],
        )
        assert set(called) == {"run_generate"}

    def test_generate_multinode_alias(self, monkeypatch):
        called = self._run(
            monkeypatch, ["orthofinder", "generate", "--multinode", "--config", "c"]
        )
        assert set(called) == {"run_generate_multinode"}

    def test_run_defaults_to_multinode(self, monkeypatch):
        called = self._run(monkeypatch, ["orthofinder", "run", "-y", "--config", "c"])
        assert set(called) == {"run_submit_multinode"}

    def test_run_single_node(self, monkeypatch):
        called = self._run(
            monkeypatch, ["orthofinder", "run", "--single-node", "-y", "--config", "c"]
        )
        assert set(called) == {"run_submit"}

    def test_conflicting_flags_exit(self, monkeypatch):
        with pytest.raises(SystemExit):
            self._run(
                monkeypatch,
                ["orthofinder", "generate", "--multinode", "--single-node"],
            )
