"""
convgeno.validation.trees — Validate phylogenetic tree inputs.

Functions to implement:
    validate_tree(path) -> list[str]
        Check that a Newick file is parseable and rooted.
        Returns a list of error messages (empty = valid).

    check_tips_match_species(tree_path, species_list) -> tuple[set, set]
        Compare tree tip labels against a species list.
        Returns (in_tree_not_in_list, in_list_not_in_tree).
        Both sets should be empty for a valid run.

    is_ultrametric(tree_path, tolerance=1e-3) -> bool
        Check whether the tree is ultrametric (required by CAFE-5).
"""
