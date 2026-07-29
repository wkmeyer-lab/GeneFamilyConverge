"""Tests for convgeno.external.r8s."""

from __future__ import annotations

import subprocess
from io import StringIO
from pathlib import Path

import pytest
from Bio import Phylo

import convgeno.utils.command_runner as command_runner
from convgeno.external import r8s
from convgeno.external.r8s import Calibration

# Simulated r8s stdout: the input tree is echoed near the top, an ASCII
# chronogram follows, and the *last* tree line (from describe
# plot=tree_description) carries the dated/ultrametric tree.
R8S_OUTPUT = (
    "[r8s banner]\n"
    "tree nj_tree = [&R] ((human:0.10,cat:0.12):0.05,dog:0.20);\n"
    "[chronogram plot ...]\n"
    "tree nj_tree = [&R] ((human:47,cat:47):47,dog:94);\n"
)

DATED_TREE = "((human:47,cat:47):47,dog:94);"


# ===================================================================
#  count_alignment_sites
# ===================================================================


class TestCountAlignmentSites:
    def test_returns_column_count(self, species_tree_alignment_file: Path):
        assert r8s.count_alignment_sites(species_tree_alignment_file) == 12

    def test_missing_file_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="alignment not found"):
            r8s.count_alignment_sites(tmp_path / "nope.fa")

    def test_ragged_alignment_raises(self, tmp_path: Path):
        p = tmp_path / "ragged.fa"
        p.write_text(">a\nMKTL\n>b\nMK\n", encoding="utf-8")
        with pytest.raises(ValueError, match="not rectangular"):
            r8s.count_alignment_sites(p)

    def test_empty_file_raises(self, tmp_path: Path):
        p = tmp_path / "empty.fa"
        p.write_text("", encoding="utf-8")
        with pytest.raises(ValueError, match="No sequences"):
            r8s.count_alignment_sites(p)


# ===================================================================
#  find_results_dir / find_species_tree_files
# ===================================================================


class TestFindSpeciesTreeFiles:
    def test_finds_tree_and_alignment(self, orthofinder_output_dir: Path):
        tree, aln = r8s.find_species_tree_files(orthofinder_output_dir)
        assert tree.name == "SpeciesTree_rooted.txt"
        assert aln.name == "SpeciesTreeAlignment.fa"
        assert tree.is_file() and aln.is_file()

    def test_accepts_results_dir_directly(self, orthofinder_output_dir: Path):
        results = orthofinder_output_dir / "Results_test"
        assert r8s.find_results_dir(results) == results

    def test_not_a_directory_raises(self, tmp_path: Path):
        with pytest.raises(NotADirectoryError):
            r8s.find_results_dir(tmp_path / "does_not_exist")

    def test_no_results_dir_raises(self, tmp_path: Path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(FileNotFoundError, match="No OrthoFinder Results"):
            r8s.find_results_dir(empty)

    def test_missing_alignment_raises_with_dendroblast_hint(self, tmp_path: Path):
        # A dendroblast-style run: species tree present, no alignment.
        root = tmp_path / "of"
        st = root / "Results_x" / "Species_Tree"
        st.mkdir(parents=True)
        (st / "SpeciesTree_rooted.txt").write_text("(a,b);", encoding="utf-8")
        with pytest.raises(FileNotFoundError, match="dendroblast"):
            r8s.find_species_tree_files(root)

    def test_picks_newest_results_dir(self, tmp_path: Path):
        import os

        root = tmp_path / "of"
        for name, mtime in (("Results_old", 1_000_000), ("Results_new", 2_000_000)):
            st = root / name / "Species_Tree"
            st.mkdir(parents=True)
            (st / "SpeciesTree_rooted.txt").write_text("(a,b);", encoding="utf-8")
            os.utime(root / name, (mtime, mtime))
        assert r8s.find_results_dir(root).name == "Results_new"

    def _make_multinode_layout(self, root: Path) -> Path:
        """Build the deep multi-node layout; return the WorkingDirectory."""
        work_dir = root / "Results_prep" / "WorkingDirectory"
        final_st = work_dir / "OrthoFinder" / "Results_final" / "Species_Tree"
        final_st.mkdir(parents=True)
        (final_st / "SpeciesTree_rooted.txt").write_text("(a,b);", encoding="utf-8")
        return work_dir

    def test_multinode_via_working_dir_pointer(self, tmp_path: Path):
        # Results are deep under the WorkingDirectory; a sibling pointer file
        # names it (as the multi-node prepare job writes).
        root = tmp_path / "orthofinder_multinode_x"
        work_dir = self._make_multinode_layout(root)
        pointer = root.parent / f"{root.name}_working_dir_path.txt"
        pointer.write_text(str(work_dir), encoding="utf-8")
        assert r8s.find_results_dir(root).name == "Results_final"

    def test_multinode_deep_recursive_without_pointer(self, tmp_path: Path):
        # Even with no pointer, the recursive fallback finds the deep Results.
        root = tmp_path / "of"
        self._make_multinode_layout(root)
        assert r8s.find_results_dir(root).name == "Results_final"


class TestRootSpanningTaxa:
    def test_returns_two_tips_spanning_root(self):
        # root children are (human,cat) and dog -> one tip from each side
        left, right = r8s.root_spanning_taxa("((human:1,cat:1):1,dog:2);")
        assert left in {"human", "cat"}
        assert right == "dog"
        assert left != right


# ===================================================================
#  sanitize_tree_for_r8s
# ===================================================================


class TestSanitizeTreeForR8s:
    def _parse(self, newick: str):
        return Phylo.read(StringIO(newick), "newick")

    def test_strips_internal_support_keeps_tips(self, species_tree_file: Path):
        clean = r8s.sanitize_tree_for_r8s(species_tree_file)
        tree = self._parse(clean)
        # No internal node retains a support/confidence value.
        assert all(nt.confidence is None for nt in tree.get_nonterminals())
        # Tips and their branch lengths survive.
        tips = tree.get_terminals()
        assert sorted(t.name for t in tips) == ["cat", "dog", "human"]
        assert all(t.branch_length is not None for t in tips)

    def test_accepts_string_and_path(self, species_tree_file: Path):
        from_path = r8s.sanitize_tree_for_r8s(species_tree_file)
        from_str = r8s.sanitize_tree_for_r8s(species_tree_file.read_text())
        for out in (from_path, from_str):
            assert out.endswith(";")
            assert "human" in out and "cat" in out and "dog" in out


# ===================================================================
#  build_r8s_control
# ===================================================================


class TestBuildR8sControl:
    TREE = "((human:0.1,cat:0.12):0.05,dog:0.2);"

    def test_point_calibration_emits_fixage(self):
        ctl = r8s.build_r8s_control(
            self.TREE, 12, [Calibration("hc", ("human", "cat"), age=94)]
        )
        assert "blformat nsites=12 lengths=persite ultrametric=no;" in ctl
        assert "mrca hc human cat;" in ctl
        assert "fixage taxon=hc age=94;" in ctl
        assert "tree nj_tree = [&R] " in ctl
        assert "divtime method=pl algorithm=tn" in ctl
        assert "describe plot=tree_description;" in ctl

    def test_default_is_single_pl_fit_no_crossv(self):
        ctl = r8s.build_r8s_control(
            self.TREE, 12, [Calibration("hc", ("human", "cat"), age=94)]
        )
        assert "set smoothing=100;" in ctl
        assert "divtime method=pl algorithm=tn;" in ctl
        assert "crossv" not in ctl  # cross-validation is opt-in

    def test_custom_smoothing(self):
        ctl = r8s.build_r8s_control(
            self.TREE, 12, [Calibration("hc", ("human", "cat"), age=94)], smoothing=5
        )
        assert "set smoothing=5;" in ctl

    def test_cross_validate_opt_in(self):
        ctl = r8s.build_r8s_control(
            self.TREE,
            12,
            [Calibration("hc", ("human", "cat"), age=94)],
            cross_validate=True,
        )
        assert "crossv=yes" in ctl
        assert "cvNum=8" in ctl
        assert "set smoothing" not in ctl  # CV selects smoothing itself

    def test_window_calibration_emits_constrain(self):
        ctl = r8s.build_r8s_control(
            self.TREE,
            12,
            [Calibration("root", ("human", "dog"), min_age=100, max_age=120)],
        )
        assert "constrain taxon=root min_age=100 max_age=120;" in ctl

    def test_multiple_calibrations(self):
        cals = [
            Calibration("hc", ("human", "cat"), age=94),
            Calibration("root", ("human", "dog"), min_age=100, max_age=120),
        ]
        ctl = r8s.build_r8s_control(self.TREE, 12, cals)
        assert "mrca hc human cat;" in ctl
        assert "mrca root human dog;" in ctl
        assert "fixage taxon=hc age=94;" in ctl
        assert "constrain taxon=root min_age=100 max_age=120;" in ctl

    def test_appends_semicolon_when_missing(self):
        ctl = r8s.build_r8s_control(
            "((a:1,b:1):1,c:2)", 10, [Calibration("ab", ("a", "b"), age=1)]
        )
        assert "tree nj_tree = [&R] ((a:1,b:1):1,c:2);" in ctl

    def test_empty_calibrations_raises(self):
        with pytest.raises(ValueError, match="At least one calibration"):
            r8s.build_r8s_control(self.TREE, 12, [])

    def test_nonpositive_nsites_raises(self):
        with pytest.raises(ValueError, match="positive integer"):
            r8s.build_r8s_control(
                self.TREE, 0, [Calibration("hc", ("human", "cat"), age=1)]
            )

    def test_calibration_without_age_raises(self):
        with pytest.raises(ValueError, match="needs a fixed"):
            r8s.build_r8s_control(self.TREE, 12, [Calibration("hc", ("human", "cat"))])

    def test_calibration_wrong_taxa_count_raises(self):
        bad = Calibration("hc", ("human",), age=1)  # only one taxon
        with pytest.raises(ValueError, match="exactly two"):
            r8s.build_r8s_control(self.TREE, 12, [bad])


# ===================================================================
#  build_command
# ===================================================================


class TestBuildCommand:
    def test_default_tool(self):
        assert r8s.build_command("ctl.txt") == ["r8s", "-b", "-f", "ctl.txt"]

    def test_custom_tool_path(self):
        assert r8s.build_command("/x/ctl", tool_path="/opt/r8s") == [
            "/opt/r8s",
            "-b",
            "-f",
            "/x/ctl",
        ]


# ===================================================================
#  parse_ultrametric_tree
# ===================================================================


class TestParseUltrametricTree:
    def test_extracts_last_tree_and_strips_rooting(self):
        assert r8s.parse_ultrametric_tree(R8S_OUTPUT) == DATED_TREE

    def test_accepts_path(self, tmp_path: Path):
        p = tmp_path / "r8s_tmp.txt"
        p.write_text(R8S_OUTPUT, encoding="utf-8")
        assert r8s.parse_ultrametric_tree(p) == DATED_TREE

    def test_no_tree_line_raises(self):
        with pytest.raises(ValueError, match="Could not find a tree"):
            r8s.parse_ultrametric_tree("some log output\nwith no tree line\n")


# ===================================================================
#  make_ultrametric (r8s mocked via command_runner)
# ===================================================================


def _fake_run_factory(captured: dict):
    def fake_run(cmd, log_dir=None, dry_run=False, **kwargs):
        captured["cmd"] = cmd
        captured["ctl"] = Path(cmd[-1]).read_text(encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, stdout=R8S_OUTPUT, stderr="")

    return fake_run


class TestMakeUltrametric:
    def test_end_to_end_writes_tree_and_stats(
        self, orthofinder_output_dir: Path, tmp_path: Path, monkeypatch
    ):
        captured: dict = {}
        monkeypatch.setattr(command_runner, "run", _fake_run_factory(captured))
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: True)

        out_tree = tmp_path / "cafe_input" / "ultrametric.nwk"
        work = tmp_path / "work"
        stats = r8s.make_ultrametric(
            orthofinder_output_dir,
            [Calibration("hc", ("human", "cat"), age=94)],
            out_tree=out_tree,
            work_dir=work,
        )

        assert out_tree.read_text(encoding="utf-8").strip() == DATED_TREE
        assert stats["nsites"] == 12
        assert stats["n_calibrations"] == 1
        assert sorted(stats["tips"]) == ["cat", "dog", "human"]
        assert (work / "r8s_ctl_file.txt").is_file()
        assert (work / "r8s_tmp.txt").is_file()
        assert captured["cmd"][:2] == ["r8s", "-b"]
        assert "blformat nsites=12 lengths=persite" in captured["ctl"]
        assert "fixage taxon=hc age=94;" in captured["ctl"]
        assert "set smoothing=100;" in captured["ctl"]  # single PL fit by default

    def test_no_calibration_uses_root_anchor(
        self, orthofinder_output_dir: Path, tmp_path: Path, monkeypatch
    ):
        # "forget about lambda": empty calibrations -> relative-time tree anchored
        # at the root, so the auto-run never blocks on a species-pair calibration.
        captured: dict = {}
        monkeypatch.setattr(command_runner, "run", _fake_run_factory(captured))
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: True)

        stats = r8s.make_ultrametric(
            orthofinder_output_dir,
            [],  # no calibrations
            out_tree=tmp_path / "t.nwk",
            work_dir=tmp_path / "w",
            root_age=5,
        )
        assert stats["relative_time"] is True
        assert stats["n_calibrations"] == 1
        # root anchored via two root-spanning tips at the given age
        assert "mrca root " in captured["ctl"]
        assert "fixage taxon=root age=5;" in captured["ctl"]

    def test_nsites_override_used_in_control(
        self, orthofinder_output_dir: Path, tmp_path: Path, monkeypatch
    ):
        captured: dict = {}
        monkeypatch.setattr(command_runner, "run", _fake_run_factory(captured))
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: True)

        r8s.make_ultrametric(
            orthofinder_output_dir,
            [Calibration("hc", ("human", "cat"), age=94)],
            out_tree=tmp_path / "t.nwk",
            work_dir=tmp_path / "w",
            nsites=999,
        )
        assert "blformat nsites=999 lengths=persite" in captured["ctl"]

    def test_dry_run_does_not_write_tree_or_need_r8s(
        self, orthofinder_output_dir: Path, tmp_path: Path, monkeypatch
    ):
        # r8s deliberately reported as unavailable; dry run must not check it.
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: False)

        out_tree = tmp_path / "dry.nwk"
        stats = r8s.make_ultrametric(
            orthofinder_output_dir,
            [Calibration("hc", ("human", "cat"), age=94)],
            out_tree=out_tree,
            work_dir=tmp_path / "w",
            dry_run=True,
        )
        assert stats["dry_run"] is True
        assert not out_tree.exists()
        assert (tmp_path / "w" / "r8s_ctl_file.txt").is_file()  # control still written

    def test_unknown_calibration_taxon_raises(
        self, orthofinder_output_dir: Path, tmp_path: Path, monkeypatch
    ):
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: True)
        monkeypatch.setattr(command_runner, "run", _fake_run_factory({}))
        with pytest.raises(ValueError, match="frog"):
            r8s.make_ultrametric(
                orthofinder_output_dir,
                [Calibration("x", ("human", "frog"), age=10)],
                out_tree=tmp_path / "t.nwk",
                work_dir=tmp_path / "w",
            )

    def test_r8s_unavailable_raises(
        self, orthofinder_output_dir: Path, tmp_path: Path, monkeypatch
    ):
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: False)
        with pytest.raises(FileNotFoundError, match="r8s"):
            r8s.make_ultrametric(
                orthofinder_output_dir,
                [Calibration("hc", ("human", "cat"), age=94)],
                out_tree=tmp_path / "t.nwk",
                work_dir=tmp_path / "w",
            )

    def test_input_tree_dates_user_tree_nsites_from_of(
        self, orthofinder_output_dir: Path, tmp_path: Path, monkeypatch
    ):
        # A user tree is dated, but nsites still comes from OrthoFinder's alignment.
        captured: dict = {}
        monkeypatch.setattr(command_runner, "run", _fake_run_factory(captured))
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: True)

        user_tree = tmp_path / "user.nwk"
        user_tree.write_text("((human:0.2,cat:0.2):0.1,dog:0.4);\n", encoding="utf-8")

        stats = r8s.make_ultrametric(
            orthofinder_output_dir,
            [Calibration("hc", ("human", "cat"), age=94)],
            out_tree=tmp_path / "t.nwk",
            work_dir=tmp_path / "w",
            input_tree=user_tree,
        )
        assert stats["nsites"] == 12  # from the fixture's 12-column alignment
        assert "0.20000000" in captured["ctl"]  # the USER tree's branch lengths

    def test_input_tree_with_nsites_needs_no_orthofinder(
        self, tmp_path: Path, monkeypatch
    ):
        # input_tree + explicit nsites: no OrthoFinder output is read at all.
        captured: dict = {}
        monkeypatch.setattr(command_runner, "run", _fake_run_factory(captured))
        monkeypatch.setattr(command_runner, "check_tool_available", lambda _t: True)

        user_tree = tmp_path / "user.nwk"
        user_tree.write_text("((human:0.2,cat:0.2):0.1,dog:0.4);\n", encoding="utf-8")

        stats = r8s.make_ultrametric(
            None,  # no OrthoFinder output dir
            [Calibration("hc", ("human", "cat"), age=94)],
            out_tree=tmp_path / "t.nwk",
            work_dir=tmp_path / "w",
            input_tree=user_tree,
            nsites=321,
        )
        assert stats["nsites"] == 321
        assert "blformat nsites=321 lengths=persite" in captured["ctl"]

    def test_input_tree_without_nsites_or_of_raises(self, tmp_path: Path):
        user_tree = tmp_path / "user.nwk"
        user_tree.write_text("((human:0.2,cat:0.2):0.1,dog:0.4);\n", encoding="utf-8")
        with pytest.raises(ValueError, match="nsites is required"):
            r8s.make_ultrametric(
                None,
                [Calibration("hc", ("human", "cat"), age=94)],
                out_tree=tmp_path / "t.nwk",
                work_dir=tmp_path / "w",
                input_tree=user_tree,
            )

    def test_missing_input_tree_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="Input species tree"):
            r8s.make_ultrametric(
                None,
                [Calibration("hc", ("human", "cat"), age=94)],
                out_tree=tmp_path / "t.nwk",
                work_dir=tmp_path / "w",
                input_tree=tmp_path / "does_not_exist.nwk",
                nsites=100,
            )
