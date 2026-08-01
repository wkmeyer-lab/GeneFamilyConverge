# make_categorical_phenotype_tree.R
#
# Build the CATEGORICAL PHENOTYPE TREE = the CAFE-5 "-y" (multi-lambda) tree.
# It has the SAME topology as the ultrametric species tree, but every branch is
# labelled with an integer naming its phenotype / rate class: tips get their
# category from the phenotype table, internal branches from ancestral-state
# reconstruction. Written to Newick with those integers in the branch-length
# slot, e.g. ((cat:2,horse:2):2,(rat:1,mouse:1):1);
#
# Because it is built ON the CAFE -t tree, the -y tree is tip-for-tip and
# node-for-node identical to -t (differing only in the branch-length slot).
#
# Reads:
#   <interim_dir>/cafe_input/species_tree.nwk   (the CAFE -t tree written by
#                                                 prepare_cafe_inputs.py;
#                                                 fallback: ultrametric.output)
#   inputs.phenotype_table                       (tip-label column + phenotype column)
#
# Writes:
#   <interim_dir>/cafe_input/lambda_tree.nwk     (the CAFE -y tree)
#   <interim_dir>/cafe_input/lambda_legend.tsv   (integer <-> category key + tip counts)
#   <interim_dir>/cafe_input/species_tree_pruned.nwk   (ONLY when unclassified: drop)
#
# Sources reusable R functions from Src/Reu/r/categorical_phenotype_tree.R
#
# Requires: ape + castor + expm (the vendored ASR engine in Src/Reu/r/anc_recon.R;
# no RERconverge needed) and yaml (only when --config is used). Optionally set
# --rer-functions / phenotype_tree.rer_functions to the lab's RERConvergeFunctions.R
# to override the vendored getAncLiks.
#
# Usage:
#   Rscript Src/Loc/scripts/make_categorical_phenotype_tree.R \
#       --config Src/Loc/configs/example_config.yaml
#   # any --flag overrides the matching config value, e.g.:
#   Rscript Src/Loc/scripts/make_categorical_phenotype_tree.R --config … \
#       --model SYM --pheno-col DietClass --unclassified drop

# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Locate this script (run via Rscript) so we can source the Src/Reu/r helpers by
# a repo-relative path rather than a hard-coded one.
script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = FALSE))
  NA_character_
}

# Minimal "--flag value" / "--flag=value" / "--flag" (boolean) parser. Flag names
# are normalised so --pheno-col and --pheno_col are equivalent.
parse_cli <- function(argv) {
  opts <- list()
  i <- 1L
  while (i <= length(argv)) {
    tok <- argv[[i]]
    if (!startsWith(tok, "--")) {
      stop(sprintf("Unexpected argument '%s' (expected --flag).", tok), call. = FALSE)
    }
    key <- sub("^--", "", tok)
    if (grepl("=", key, fixed = TRUE)) {
      pos <- regexpr("=", key, fixed = TRUE)
      name <- substr(key, 1, pos - 1)
      opts[[gsub("-", "_", name)]] <- substr(key, pos + 1, nchar(key))
      i <- i + 1L
    } else if (i + 1L <= length(argv) && !startsWith(argv[[i + 1L]], "--")) {
      opts[[gsub("-", "_", key)]] <- argv[[i + 1L]]
      i <- i + 2L
    } else {
      opts[[gsub("-", "_", key)]] <- TRUE
      i <- i + 1L
    }
  }
  opts
}

load_yaml <- function(path) {
  if (is.null(path) || isTRUE(is.na(path)) || !nzchar(path)) return(list())
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' package is required to read --config; install it or pass ",
         "all values as --flags.", call. = FALSE)
  }
  yaml::read_yaml(path)
}

# Nested lookup: dig(cfg, "a", "b") -> cfg$a$b, or NULL if any level is absent.
dig <- function(x, ...) {
  for (k in c(...)) {
    if (!is.list(x) || is.null(x[[k]])) return(NULL)
    x <- x[[k]]
  }
  x
}

# YAML list-of-pairs [[from,to],[a,b]] -> list(c("from","to"), c("a","b")); [] -> NULL.
norm_pairs <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  lapply(x, function(p) as.character(unlist(p)))
}

is_unset <- function(x) {
  is.null(x) || isTRUE(is.na(x)) ||
    (is.character(x) && length(x) == 1 && x %in% c("", "null", "NULL", "none", "None"))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

run <- function(argv) {
  opts <- parse_cli(argv)
  cfg <- load_yaml(opts$config)

  # Respect phenotype_tree.enabled so the workflows can call this step
  # unconditionally (like the other pipeline steps). Absent/true -> run.
  enabled <- opts$enabled %||% dig(cfg, "phenotype_tree", "enabled")
  if (isFALSE(enabled) ||
      (is.character(enabled) && tolower(enabled) %in% c("false", "no", "0", "off"))) {
    cat("phenotype_tree.enabled is false; skipping the categorical phenotype tree step.\n")
    return(invisible(0))
  }

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("The 'ape' package is required.", call. = FALSE)
  }

  # --- source the reusable helpers -----------------------------------------
  reu_dir <- opts$reu_dir
  if (is.null(reu_dir)) {
    sp <- script_path()
    if (is.na(sp)) {
      stop("Could not locate this script; pass --reu-dir <Src/Reu/r>.", call. = FALSE)
    }
    # this script is Src/Loc/scripts/<file>.R -> repo root is three levels up
    repo_root <- normalizePath(file.path(dirname(sp), "..", "..", ".."), mustWork = FALSE)
    reu_dir <- file.path(repo_root, "Src", "Reu", "r")
  }
  helper <- file.path(reu_dir, "categorical_phenotype_tree.R")
  if (!file.exists(helper)) {
    stop(sprintf("Cannot find helper %s (set --reu-dir).", helper), call. = FALSE)
  }
  source(helper)

  # Vendored ancestral-state-reconstruction engine (getAncLiks); makes the step
  # depend on castor + expm instead of the whole RERconverge package. The optional
  # --rer-functions below can override this with the lab's RERConvergeFunctions.R.
  anc <- file.path(reu_dir, "anc_recon.R")
  if (!file.exists(anc)) {
    stop(sprintf("Cannot find %s (the vendored ASR engine; set --reu-dir).", anc),
         call. = FALSE)
  }
  source(anc)

  # Optionally source the lab's RERConvergeFunctions.R so getAncLiks() is
  # available when the RERconverge package is not installed.
  rer_functions <- opts$rer_functions %||% dig(cfg, "phenotype_tree", "rer_functions")
  if (!is_unset(rer_functions)) {
    if (!file.exists(rer_functions)) {
      stop(sprintf("--rer-functions file not found: %s", rer_functions), call. = FALSE)
    }
    source(rer_functions)
  }

  # --- resolve inputs (CLI overrides config overrides default) -------------
  interim <- dig(cfg, "outputs", "interim_dir") %||% "Data/interim"
  cafe_dir <- file.path(interim, "cafe_input")

  tree_path <- opts$tree
  if (is.null(tree_path)) {
    canonical <- file.path(cafe_dir, "species_tree.nwk")
    tree_path <- if (file.exists(canonical)) canonical else dig(cfg, "ultrametric", "output")
  }
  if (is.null(tree_path)) {
    stop("No species tree: pass --tree, or run prepare_cafe_inputs first, or set ",
         "ultrametric.output in --config.", call. = FALSE)
  }
  if (!file.exists(tree_path)) stop(sprintf("Tree not found: %s", tree_path), call. = FALSE)

  pheno_path <- opts$phenotypes %||% dig(cfg, "inputs", "phenotype_table")
  if (is.null(pheno_path)) {
    stop("No phenotype table: pass --phenotypes or set inputs.phenotype_table.", call. = FALSE)
  }
  if (!file.exists(pheno_path)) {
    stop(sprintf("Phenotype table not found: %s", pheno_path), call. = FALSE)
  }

  id_col     <- opts$id_col     %||% dig(cfg, "phenotype_tree", "id_col")     %||% "species"
  pheno_col  <- opts$pheno_col  %||% dig(cfg, "phenotype_tree", "pheno_col")  %||% "phenotype"
  model      <- opts$model      %||% dig(cfg, "phenotype_tree", "model")      %||% "ER"
  root_prior <- opts$root_prior %||% dig(cfg, "phenotype_tree", "root_prior") %||% "auto"

  anctrait <- opts$anctrait %||% dig(cfg, "phenotype_tree", "anctrait")
  if (is_unset(anctrait)) anctrait <- NULL

  unclassified <- opts$unclassified %||% dig(cfg, "phenotype_tree", "unclassified") %||%
    "background"
  if (!unclassified %in% c("background", "drop")) {
    stop(sprintf("unclassified must be 'background' or 'drop' (got '%s').", unclassified),
         call. = FALSE)
  }
  background_label <- opts$background_label %||%
    dig(cfg, "phenotype_tree", "background_label") %||% "background"

  merges <- dig(cfg, "phenotype_tree", "merges")
  substitutions <- norm_pairs(dig(merges, "substitutions"))
  merge_only    <- norm_pairs(dig(merges, "merge_only"))

  out_tree <- opts$out    %||% dig(cfg, "phenotype_tree", "output") %||%
    file.path(cafe_dir, "lambda_tree.nwk")
  out_legend <- opts$legend %||% dig(cfg, "phenotype_tree", "legend") %||%
    file.path(cafe_dir, "lambda_legend.tsv")

  # --- read the -t tree ----------------------------------------------------
  tree <- ape::read.tree(tree_path)
  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop(sprintf("Could not read a single Newick tree from %s", tree_path), call. = FALSE)
  }
  if (!ape::is.rooted(tree)) {
    stop("Species tree is not rooted; CAFE-5 needs a rooted -t/-y tree.", call. = FALSE)
  }
  if (!ape::is.binary(tree)) {
    stop("Species tree is not strictly binary; CAFE-5 needs a binary tree.", call. = FALSE)
  }

  # --- read the phenotype table -------------------------------------------
  df <- utils::read.delim(pheno_path, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE,
                          na.strings = c("NA", ""), quote = "")
  for (col in c(id_col, pheno_col)) {
    if (!col %in% names(df)) {
      stop(sprintf("Column '%s' not in %s. Available: %s",
                   col, pheno_path, paste(names(df), collapse = ", ")), call. = FALSE)
    }
  }
  ids <- trimws(as.character(df[[id_col]]))
  vals <- trimws(as.character(df[[pheno_col]]))
  keep <- !is.na(ids) & nzchar(ids)
  ids <- ids[keep]; vals <- vals[keep]

  # Duplicate ids are fine if they agree; conflicting labels are an error.
  if (anyDuplicated(ids)) {
    by_id <- split(vals, ids)
    conflict <- names(Filter(function(v) length(unique(v[!is.na(v)])) > 1, by_id))
    if (length(conflict)) {
      stop(sprintf("Species with conflicting phenotype labels: %s",
                   paste(utils::head(conflict, 10), collapse = ", ")), call. = FALSE)
    }
    vals <- vapply(by_id, function(v) v[!is.na(v)][1] %||% NA_character_, character(1))
    ids <- names(by_id)
  }
  pheno <- stats::setNames(vals, ids)

  # Collapse fine-grained labels into the analysis category set.
  pheno <- merge_phenotype_categories(pheno, substitutions = substitutions,
                                      merge_only = merge_only)

  # --- assemble a complete per-tip vector, applying the unclassified policy -
  work_tree <- tree
  tip_pheno <- stats::setNames(unname(pheno[work_tree$tip.label]), work_tree$tip.label)
  missing <- names(tip_pheno)[is.na(tip_pheno)]
  extra <- setdiff(names(pheno), tree$tip.label)

  if (length(missing)) {
    if (unclassified == "background") {
      tip_pheno[missing] <- background_label
      message(sprintf("%d tip(s) without a phenotype assigned to '%s'.",
                      length(missing), background_label))
    } else {  # drop
      work_tree <- ape::drop.tip(tree, missing)
      tip_pheno <- stats::setNames(unname(pheno[work_tree$tip.label]), work_tree$tip.label)
      pruned_t <- file.path(dirname(out_tree), "species_tree_pruned.nwk")
      dir.create(dirname(pruned_t), recursive = TRUE, showWarnings = FALSE)
      ape::write.tree(work_tree, pruned_t)
      warning(sprintf(
        paste0("Dropped %d unclassified tip(s). Wrote a matching pruned -t to %s. ",
               "You MUST subset the CAFE count matrix (cafe_input.tsv) to the same ",
               "%d tips before running CAFE."),
        length(missing), pruned_t, length(work_tree$tip.label)), call. = FALSE)
    }
  }
  if (length(extra)) {
    message(sprintf("%d phenotype-table row(s) not in the tree were ignored.", length(extra)))
  }

  # --- build the categorical phenotype tree --------------------------------
  built <- build_categorical_phenotype_tree(
    tree = work_tree,
    phenotypes = tip_pheno,
    model = model,
    root_prior = root_prior,
    anctrait = anctrait
  )
  ytree <- built$tree
  legend <- built$legend
  k <- built$n_states

  # --- validate before writing --------------------------------------------
  if (!setequal(ytree$tip.label, work_tree$tip.label)) {
    stop("Internal error: -y tips differ from -t tips.", call. = FALSE)
  }
  all_equal_phylo <- get("all.equal.phylo", envir = asNamespace("ape"))
  if (!isTRUE(all_equal_phylo(work_tree, ytree, use.edge.length = FALSE))) {
    stop("Internal error: -y topology differs from -t topology.", call. = FALSE)
  }
  uniq_ints <- sort(unique(as.integer(round(ytree$edge.length))))
  if (!identical(uniq_ints, seq_len(k))) {
    warning(sprintf(
      paste0("Branch classes present {%s} != 1..%d; a category may be absent from ",
             "all branches, leaving CAFE with an unused lambda."),
      paste(uniq_ints, collapse = ", "), k), call. = FALSE)
  }

  # --- write outputs -------------------------------------------------------
  dir.create(dirname(out_tree), recursive = TRUE, showWarnings = FALSE)
  write_categorical_phenotype_tree(ytree, out_tree)

  counts <- table(tip_pheno)
  legend$n_tips <- as.integer(counts[legend$category])
  legend$n_tips[is.na(legend$n_tips)] <- 0L
  dir.create(dirname(out_legend), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(legend, out_legend, sep = "\t", quote = FALSE, row.names = FALSE)

  # --- summary -------------------------------------------------------------
  anc_note <- if (is.null(anctrait)) "" else sprintf(", anctrait=%s", anctrait)
  cat(sprintf("tree (-t):        %s\n", tree_path))
  cat(sprintf("phenotypes:       %s  [%s -> %s]\n", pheno_path, id_col, pheno_col))
  cat(sprintf("model:            %s (root_prior=%s)%s\n", model, root_prior, anc_note))
  cat(sprintf("tips:             %d\n", length(ytree$tip.label)))
  cat(sprintf("phenotype classes (%d):\n", k))
  for (r in seq_len(nrow(legend))) {
    cat(sprintf("  %d = %-24s (%d tips)\n", legend$integer[r], legend$category[r],
                legend$n_tips[r]))
  }
  cat(sprintf("categorical phenotype tree (-y) -> %s\n", out_tree))
  cat(sprintf("legend -> %s\n", out_legend))
  invisible(0)
}

if (sys.nframe() == 0) {
  tryCatch(
    run(commandArgs(trailingOnly = TRUE)),
    error = function(e) {
      message("ERROR: ", conditionMessage(e))
      quit(status = 1, save = "no")
    }
  )
}
