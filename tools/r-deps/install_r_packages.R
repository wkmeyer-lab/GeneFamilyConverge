#!/usr/bin/env Rscript
# install_r_packages.R
#
# Install the small R stack the CATEGORICAL PHENOTYPE TREE step needs:
#   Src/Loc/scripts/make_categorical_phenotype_tree.R  ->  the CAFE-5 "-y" tree.
#
# Packages (all on CRAN, no root required — installs to your user library):
#   ape    : Newick/phylo I/O + reorder() used by the ancestral reconstruction.
#   castor : fit_mk() — the Mk-model fit at the heart of the reconstruction.
#   expm   : expm()   — matrix exponential used in the pruning passes.
#   yaml   : read the pipeline --config file.
# (Matrix ships with R; phangorn and RERconverge are NOT needed — the ASR engine
#  getAncLiks() is vendored in Src/Reu/r/anc_recon.R with castor + expm only.)
#
# The repo keeps R + packages OUT of the conda env (like tools/r8s), so R is
# supplied by the system / an HPC module. Usage:
#
#   module load R            # version is cluster-specific; skip on a workstation
#   Rscript tools/r-deps/install_r_packages.R
#
# Optional overrides:
#   Rscript tools/r-deps/install_r_packages.R /path/to/lib   # target library
#   R_CRAN_MIRROR=https://cran.r-project.org Rscript tools/r-deps/install_r_packages.R
#
# conda-forge alternative (adds R to an env — keep it separate from convgeno):
#   conda install -c conda-forge r-base r-ape r-castor r-expm r-yaml

required <- c("ape", "castor", "expm", "yaml")

mirror <- Sys.getenv("R_CRAN_MIRROR", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = mirror))

# --- choose a writable install library -----------------------------------
args <- commandArgs(trailingOnly = TRUE)
explicit_lib <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else NA_character_

pick_library <- function(explicit) {
  if (!is.na(explicit)) {
    dir.create(explicit, recursive = TRUE, showWarnings = FALSE)
    return(explicit)
  }
  first <- .libPaths()[1]
  if (file.access(first, mode = 2) == 0) return(first)   # default lib is writable
  # Otherwise set up a personal library (typical on HPC / read-only system libs).
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib) || user_lib == "NULL") {
    user_lib <- file.path(path.expand("~"), "R", "library")
  } else {
    user_lib <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]][1]
  }
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  user_lib
}

lib <- pick_library(explicit_lib)

cat(sprintf("R %s\n", getRversion()))
cat(sprintf("CRAN mirror:     %s\n", mirror))
cat(sprintf("Install library: %s\n\n", lib))

# --- install only what's missing (idempotent) ----------------------------
have <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
missing <- required[!have]
if (length(missing) == 0) {
  cat("All required packages are already installed.\n")
} else {
  cat("Installing:", paste(missing, collapse = ", "), "\n\n")
  install.packages(missing, lib = lib)
}

# --- verify --------------------------------------------------------------
cat("\nVerification:\n")
ok <- TRUE
for (p in required) {
  present <- requireNamespace(p, quietly = TRUE)
  ok <- ok && present
  ver <- if (present) as.character(utils::packageVersion(p)) else "-"
  cat(sprintf("  %-8s %-8s %s\n", p, if (present) "OK" else "MISSING", ver))
}

if (!ok) {
  cat("\nERROR: some packages failed to install (see messages above). If the\n")
  cat("install library was not writable, pass one explicitly, e.g.:\n")
  cat("  Rscript tools/r-deps/install_r_packages.R \"$HOME/R/library\"\n")
  quit(status = 1, save = "no")
}

cat("\nAll set. Build the categorical phenotype tree (CAFE -y) with:\n")
cat("  Rscript Src/Loc/scripts/make_categorical_phenotype_tree.R --config <config.yaml>\n")
