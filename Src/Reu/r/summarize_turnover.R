# summarize_turnover.R
#
# Functions for summarizing gene family turnover patterns.
#
# Functions to implement:
#
#   summarize_turnover_by_branch(cafe_changes, tree)
#     Summarize total expansions and contractions per branch.
#
#   overlap_with_phenotype_transitions(cafe_changes, phenotype_tree)
#     Compare branches with high gene family turnover against
#     branches where the phenotype transitions.
#     This is one of the four association metrics from Michael's
#     thesis proposal (Chapter 4, Figure 6, step 10).
#
#   summarize_family_dynamics(cafe_counts)
#     For each family: min/max/mean copy number, total gain, total loss.
