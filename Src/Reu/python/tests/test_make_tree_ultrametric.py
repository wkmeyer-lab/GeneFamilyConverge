"""Tests for the Src/Loc/scripts/make_tree_ultrametric.py CLI (design step 3).

The script is a standalone Loc entry point (not part of the ``convgeno`` package),
so it is loaded from its file path via importlib.
"""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import convgeno.utils.command_runner as command_runner

_SCRIPT = (
    Path(__file__).resolve().parents[4]
    / "Src"
    / "Loc"
    / "scripts"
    / "make_tree_ultrametric.py"
)
_spec = importlib.util.spec_from_file_location("make_tree_ultrametric", _SCRIPT)
mtu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mtu)

# Every root-to-tip path sums to 94 -> ultrametric.
ULTRAMETRIC = "((human:47,cat:47):47,dog:94);\n"
# Unequal root-to-tip depths -> additive / non-ultrametric.
ADDITIVE = "((human:0.10,cat:0.12):0.05,dog:0.20);\n"

R8S_OUTPUT = (
    "tree nj_tree = [&R] ((human:0.10,cat:0.12):0.05,dog:0.20);\n"
    "tree nj_tree = [&R] ((human:47,cat:47):47,dog:94);\n"
)


class TestAssumeUltrametric:
    def test_copies_ultrametric_tree_and_needs_no_orthofinder(self, tmp_path: Path):
        tree = tmp_path / "dated.nwk"
        tree.write_text(ULTRAMETRIC, encoding="utf-8")
        out = tmp_path / "cafe_input" / "species_tree_ultrametric.nwk"

        rc = mtu.main(
            ["--input-tree", str(tree), "--assume-ultrametric", "-o", str(out)]
        )

        assert rc == 0
        assert out.is_file()
        text = out.read_text(encoding="utf-8")
        assert "human" in text and "dog" in text

    def test_rejects_non_ultrametric_tree(self, tmp_path: Path):
        tree = tmp_path / "additive.nwk"
        tree.write_text(ADDITIVE, encoding="utf-8")
        out = tmp_path / "out.nwk"

        rc = mtu.main(
            ["--input-tree", str(tree), "--assume-ultrametric", "-o", str(out)]
        )

        assert rc == 1
        assert not out.exists()

    def test_requires_input_tree(self, tmp_path: Path):
        rc = mtu.main(["--assume-ultrametric", "-o", str(tmp_path / "out.nwk")])
        assert rc == 1


class TestInputTreeWithR8s:
    def test_user_tree_plus_nsites_skips_orthofinder(
        self, tmp_path: Path, monkeypatch
    ):
        captured: dict = {}

        def fake_run(cmd, log_dir=None, dry_run=False, **kwargs):
            captured["ctl"] = Path(cmd[-1]).read_text(encoding="utf-8")
            return subprocess.CompletedProcess(cmd, 0, stdout=R8S_OUTPUT, stderr="")

        monkeypatch.setattr(command_runner, "run", fake_run)
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: True)

        tree = tmp_path / "additive.nwk"
        tree.write_text(ADDITIVE, encoding="utf-8")
        out = tmp_path / "out.nwk"

        # No OrthoFinder output dir given -> success proves none was read.
        rc = mtu.main(
            [
                "--input-tree",
                str(tree),
                "--nsites",
                "321",
                "-o",
                str(out),
                "--calibration",
                "hc:human,cat:94",
            ]
        )

        assert rc == 0
        assert out.is_file()
        assert "blformat nsites=321 lengths=persite" in captured["ctl"]

    def test_input_tree_without_nsites_or_of_errors(self, tmp_path: Path):
        tree = tmp_path / "additive.nwk"
        tree.write_text(ADDITIVE, encoding="utf-8")
        out = tmp_path / "out.nwk"

        rc = mtu.main(["--input-tree", str(tree), "-o", str(out)])

        assert rc == 1
        assert not out.exists()


class TestTipMatchHelpers:
    def test_species_from_dir_reads_fasta_stems(self, tmp_path: Path):
        prot = tmp_path / "proteomes"
        prot.mkdir()
        (prot / "human.faa").write_text(">x\nMK\n", encoding="utf-8")
        (prot / "cat.fasta").write_text(">x\nMK\n", encoding="utf-8")
        (prot / "README.txt").write_text("ignore me\n", encoding="utf-8")

        assert mtu._species_from_dir(prot) == {"human", "cat"}

    def test_report_tip_match_without_dir_is_noop(self, tmp_path: Path):
        tree = tmp_path / "dated.nwk"
        tree.write_text(ULTRAMETRIC, encoding="utf-8")
        # Must not raise when no species directory is available.
        mtu._report_tip_match(tree, None)
