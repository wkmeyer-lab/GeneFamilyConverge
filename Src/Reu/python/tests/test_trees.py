"""Tests for convgeno.validation.trees."""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.validation.trees import (
    check_tips_match_species,
    is_ultrametric,
    root_to_tip_depths,
    ultrametric_deviation,
    validate_tree,
)


class TestValidateTree:
    def test_valid_binary_rooted_tree(self, ultrametric_tree_file: Path):
        assert validate_tree(ultrametric_tree_file) == []

    def test_accepts_newick_string(self):
        assert validate_tree("((a:1,b:1):1,c:2);") == []

    def test_polytomy_reported(self):
        errors = validate_tree("((a:1,b:1,c:1):1,d:2);")
        assert any("polytom" in e.lower() for e in errors)

    def test_unrooted_root_reported(self):
        # A 3-way split at the root is not a rooted bifurcating tree.
        errors = validate_tree("(a:1,b:1,c:1);")
        assert any("root" in e.lower() for e in errors)

    def test_unparseable_reported(self):
        errors = validate_tree("((a,b);")  # missing close paren
        assert errors and "parse" in errors[0].lower()


class TestIsUltrametric:
    def test_true_for_dated_tree(self, ultrametric_tree_file: Path):
        assert is_ultrametric(ultrametric_tree_file) is True

    def test_false_for_additive_tree(self, species_tree_file: Path):
        # OrthoFinder-style ML tree: unequal root-to-tip distances.
        assert is_ultrametric(species_tree_file) is False

    def test_true_via_newick_string(self):
        assert is_ultrametric("((a:47,b:47):47,c:94);") is True

    def test_tolerance_is_relative_to_height(self):
        # Tiny deviation (94 vs 94.01) is within the default 0.1% tolerance.
        assert is_ultrametric("((a:47,b:47):47,c:94.01);") is True
        # A gross deviation is not.
        assert is_ultrametric("((a:47,b:47):47,c:150);") is False


class TestRootToTipDepths:
    def test_equal_depths_for_ultrametric(self, ultrametric_tree_file: Path):
        depths = root_to_tip_depths(ultrametric_tree_file)
        assert set(depths) == {"human", "cat", "dog"}
        assert max(depths.values()) - min(depths.values()) == pytest.approx(0.0)

    def test_unequal_depths_for_additive(self, species_tree_file: Path):
        depths = root_to_tip_depths(species_tree_file)
        assert set(depths) == {"human", "cat", "dog"}
        assert max(depths.values()) - min(depths.values()) > 0


class TestUltrametricDeviation:
    def test_zero_for_ultrametric(self, ultrametric_tree_file: Path):
        assert ultrametric_deviation(ultrametric_tree_file) == pytest.approx(0.0)

    def test_positive_for_additive(self, species_tree_file: Path):
        assert ultrametric_deviation(species_tree_file) > 0

    def test_via_newick_string(self):
        assert ultrametric_deviation("((a:47,b:47):47,c:94);") == pytest.approx(0.0)
        assert ultrametric_deviation("((a:47,b:47):47,c:150);") == pytest.approx(56.0)


class TestCheckTipsMatchSpecies:
    def test_exact_match(self, ultrametric_tree_file: Path):
        only_tree, only_list = check_tips_match_species(
            ultrametric_tree_file, ["human", "cat", "dog"]
        )
        assert only_tree == set()
        assert only_list == set()

    def test_reports_mismatches(self, ultrametric_tree_file: Path):
        only_tree, only_list = check_tips_match_species(
            ultrametric_tree_file, ["human", "cat", "frog"]
        )
        assert only_tree == {"dog"}  # in tree, not in list
        assert only_list == {"frog"}  # in list, not in tree
