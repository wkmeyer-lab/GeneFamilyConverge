# read_badirate_outputs.R
#
# Functions to parse BadiRate result files into clean R data structures.
#
# Functions to implement:
#
#   read_badirate_results(path)
#     Read BadiRate output -> data.frame with columns:
#     family_id, null_aic, alt_aic, delta_aic, better_model
#
#   get_phenotype_linked_families(results_df, aic_threshold = 2)
#     Filter to gene families where the phenotype-driven model
#     fits better than the null (delta_AIC > threshold).
#     Convention: delta_AIC > 2 is "substantially better".
