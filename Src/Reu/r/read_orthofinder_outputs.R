# read_orthofinder_outputs.R
#
# Functions to parse OrthoFinder result files into clean R data structures.
#
# Functions to implement:
#
#   read_gene_counts(path)
#     Read Orthogroups_GeneCount.tsv -> data.frame
#     Columns: orthogroup_id, species1, species2, ..., total
#
#   read_orthogroups(path)
#     Read Orthogroups.tsv -> named list of orthogroup contents
#
#   read_single_copy_orthologs(path)
#     Read Orthogroups_SingleCopyOrthologues.txt -> character vector
#
#   read_species_tree(path)
#     Read Species_Tree/SpeciesTree_rooted.txt -> ape::phylo object
