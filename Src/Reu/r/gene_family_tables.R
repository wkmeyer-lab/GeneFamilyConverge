# gene_family_tables.R
#
# Functions for merging and summarizing gene family data across tools.
#
# Functions to implement:
#
#   merge_family_results(orthofinder_counts, cafe_pvalues, badirate_results)
#     Join results from all three tools into a single summary table
#     keyed by orthogroup/family ID.
#
#   add_phenotype_column(family_table, phenotype_table)
#     Attach phenotype metadata to the family summary.
#
#   compute_mean_copy_number_by_phenotype(family_table, phenotype_col)
#     For each family, compute mean copy number within each phenotype group.
