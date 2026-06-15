# read_cafe_outputs.R
#
# Functions to parse CAFE-5 result files into clean R data structures.
#
# Functions to implement:
#
#   read_cafe_counts(path)
#     Read Base_count.tab -> data.frame of gene family sizes at each node
#
#   read_cafe_changes(path)
#     Read Base_change.tab -> data.frame of size changes per branch
#
#   read_cafe_family_pvalues(path)
#     Read Gamma_family_results.txt -> data.frame with family-level p-values
#
#   get_significant_families(pvalue_df, threshold = 0.05)
#     Filter to rapidly evolving gene families
#
#   read_cafe_branch_probs(path)
#     Read Base_branch_probabilities.tab -> data.frame
