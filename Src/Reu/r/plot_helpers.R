# plot_helpers.R
#
# Reusable plotting functions (ggplot2-based).
# Loc scripts call these; this file contains no project-specific paths.
#
# Functions to implement:
#
#   plot_copy_number_heatmap(family_table, species_order)
#     Heatmap of gene family copy number across species.
#
#   plot_family_size_distribution(family_table)
#     Histogram of family sizes across all orthogroups.
#
#   plot_turnover_on_tree(tree, branch_turnover)
#     Phylogeny with branches colored by expansion/contraction.
#     Uses ape or ggtree.
#
#   plot_phenotype_comparison(family_table, phenotype_col, metric_col)
#     Boxplot comparing a metric (e.g., copy number) across phenotype groups.
