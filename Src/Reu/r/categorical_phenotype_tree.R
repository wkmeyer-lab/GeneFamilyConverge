# categorical_phenotype_tree.R
#
# Reusable, path-free helpers to build a CATEGORICAL PHENOTYPE TREE and serialise
# it as a CAFE-5 rate ("-y" / lambda) tree.
#
# A categorical phenotype tree is the species-tree topology with EVERY branch
# (tips and internal branches) labelled by an integer naming the phenotype
# category on that branch. The tip categories come from the phenotype table; the
# internal-branch categories come from ancestral-state reconstruction. Written to
# Newick with those integers in the branch-length slot, the result is exactly
# CAFE-5's multi-lambda "-y" tree, e.g.
#     ((cat:2,horse:2):2,(rat:1,mouse:1):1);
# where the integers are lambda (rate-class) ids, not branch lengths.
#
# WHY WE DO NOT CALL char2TreeCategorical() DIRECTLY
# --------------------------------------------------
# The Meyer-lab RERconverge function char2TreeCategorical() is the natural
# reuse target, and this file DOES reuse its ancestral-state-reconstruction
# engine, RERconverge::getAncLiks() (char2TreeCategorical calls exactly the same
# routine internally). But char2TreeCategorical() ends with unroot(), returning
# an UNROOTED tree (basal trifurcation). CAFE-5 requires a ROOTED, strictly
# binary "-y" tree whose topology is identical to the ultrametric "-t" species
# tree. So we run getAncLiks() on the rooted "-t" tree and assign the integer
# states onto that same rooted topology -- keeping "-t" and "-y" tip-for-tip and
# node-for-node identical (differing only in the branch-length slot).
#
# DEPENDENCIES (must be available in the R environment; none are hard-coded paths)
#   - RERconverge  : getAncLiks() (which itself pulls castor::fit_mk, phangorn, expm)
#   - ape          : read.tree / write.tree / drop.tip / phylo handling
#   - categoricalDropTip() from CategoricalDropTip.R -- ONLY for
#     prune_categorical_phenotype_tree(); the CALLER sources that file (kept out
#     of here so this file stays path-free per the Src/Reu convention).
#
# All inputs are passed as arguments (trees, vectors, paths); nothing project- or
# cluster-specific lives here.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.cpt_require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}

# Resolve getAncLiks whether RERconverge is installed-and-attached or its
# functions were source()d into the global environment.
.cpt_get_anc_liks <- function() {
  if (requireNamespace("RERconverge", quietly = TRUE)) {
    ns <- asNamespace("RERconverge")
    if (exists("getAncLiks", envir = ns, inherits = FALSE)) {
      return(get("getAncLiks", envir = ns))
    }
  }
  if (exists("getAncLiks", mode = "function")) {
    return(get("getAncLiks", mode = "function"))
  }
  stop(
    "getAncLiks() not found. Install/attach RERconverge, or source the lab's ",
    "RERConvergeFunctions.R so getAncLiks() is available.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. Collapse fine-grained phenotype labels into the analysis category set
# ---------------------------------------------------------------------------
#
# vec           : named character vector (names = tip labels, values = raw label).
# substitutions : list of c(from, to); every value exactly equal to `from`
#                 becomes `to` (a whole-category relabel/collapse, e.g.
#                 c("Piscivore", "Carnivore")). This is the '=' semantics of the
#                 'u=' pass in MakeCategoricalPhenotypeTree.R, made exact-match
#                 (not gsub) on purpose: gsub silently mangles substrings, which
#                 is unsafe for clean category names.
# merge_only    : list of c(a, b); only HYBRID values that contain BOTH tokens
#                 (e.g. "Carnivore/Piscivore") become `b`; a standalone `a` is
#                 left unchanged. Mirrors the merge-only ('o=') pass.
#
# Returns the modified vector with names preserved.
merge_phenotype_categories <- function(vec, substitutions = NULL, merge_only = NULL) {
  out <- vec
  if (!is.null(merge_only)) {
    for (pair in merge_only) {
      a <- pair[[1]]; b <- pair[[2]]
      hybrid <- grepl(a, out, fixed = TRUE) & grepl(b, out, fixed = TRUE)
      out[hybrid] <- b
    }
  }
  if (!is.null(substitutions)) {
    for (pair in substitutions) {
      from <- pair[[1]]; to <- pair[[2]]
      out[which(out == from)] <- to   # which() drops NAs (unclassified tips)
    }
  }
  out
}

# ---------------------------------------------------------------------------
# 2. Build the categorical phenotype tree on a ROOTED tree
# ---------------------------------------------------------------------------
#
# tree       : a rooted, binary ape::phylo (typically the ultrametric CAFE "-t"
#              species tree). Its (time) branch lengths are used by the ASR.
# phenotypes : named character vector; names = tip labels (must cover EVERY tip
#              of `tree`, since -y and -t must share all tips), values = category.
# model      : rate model for getAncLiks: "ER", "SYM", or "ARD".
# root_prior : root-state prior for getAncLiks ("auto", "flat", "empirical", ...).
# anctrait   : if non-NULL, force EVERY internal branch to this category instead
#              of inferring ancestral states (tip branches keep their own
#              category). Must be one of the phenotype categories. Mirrors the
#              'g=' option of the lab script.
#
# Returns list(
#   tree     = <phylo, same rooted topology as `tree`, edge.length = integer
#               category code of each branch's child node>,
#   legend   = data.frame(integer, category)  -- the code<->category key
#               (char2TreeCategorical only prints this; here it is returned),
#   n_states = number of categories (k)
# )
#
# Encoding matches char2TreeCategorical's edge.length = states[edge[,2]] scheme
# (each branch carries the integer state of the node BELOW it), but is 1..k for
# every k (no 0/1 special-case), which is what CAFE-5 expects.
build_categorical_phenotype_tree <- function(tree, phenotypes,
                                             model = "ER",
                                             root_prior = "auto",
                                             anctrait = NULL) {
  .cpt_require_pkg("ape")

  if (!inherits(tree, "phylo")) stop("`tree` must be an ape::phylo object.", call. = FALSE)
  if (is.null(names(phenotypes))) {
    stop("`phenotypes` must be a named vector (names = tip labels).", call. = FALSE)
  }

  tips <- tree$tip.label
  missing_tips <- setdiff(tips, names(phenotypes))
  if (length(missing_tips) > 0) {
    stop(sprintf(
      paste0("Every tip needs a phenotype (a CAFE -y tree must share all tips ",
             "with -t). %d tip(s) missing, e.g.: %s"),
      length(missing_tips), paste(utils::head(missing_tips, 10), collapse = ", ")),
      call. = FALSE)
  }

  # Order the phenotype vector to tip-node order (nodes 1..Ntip).
  pheno <- as.character(phenotypes[tips])
  names(pheno) <- tips
  if (any(is.na(pheno) | pheno == "")) {
    bad <- tips[is.na(pheno) | pheno == ""]
    stop(sprintf("Tip(s) with empty/NA phenotype: %s",
                 paste(utils::head(bad, 10), collapse = ", ")), call. = FALSE)
  }

  # Deterministic category -> integer key (1..k), sorted by category name so the
  # mapping is reproducible across runs.
  categories <- sort(unique(pheno))
  k <- length(categories)
  if (k < 2) {
    stop("Need at least 2 phenotype categories to build a rate tree.", call. = FALSE)
  }
  name2index <- stats::setNames(seq_along(categories), categories)
  legend <- data.frame(
    integer  = unname(name2index),
    category = names(name2index),
    stringsAsFactors = FALSE
  )

  ytree <- tree
  ytree$node.label <- NULL   # CAFE rejects internal node labels
  ytree$root.edge  <- NULL   # ... and a stray root edge
  Ntip <- length(tips)

  if (is.null(anctrait)) {
    # ---- infer ancestral states with the lab's ASR engine, on the ROOTED tree
    getAncLiks <- .cpt_get_anc_liks()
    mapped <- unname(name2index[pheno])   # integer tip states, in tip-node order
    ancliks <- getAncLiks(tree, mapped, rate_model = model, root_prior = root_prior)
    # ancliks: one row per internal node, ordered node (Ntip+1), (Ntip+2), ...;
    # columns are states 1..k in sorted order (== our codes).
    internal <- apply(ancliks, 1, which.max)
    states <- c(mapped, internal)              # indexed by node number
    ytree$edge.length <- states[ytree$edge[, 2]]  # each branch = child node state
  } else {
    # ---- force all internal branches to the ancestral class (no ASR)
    if (!anctrait %in% categories) {
      stop(sprintf("anctrait '%s' is not one of the categories: %s",
                   anctrait, paste(categories, collapse = ", ")), call. = FALSE)
    }
    j <- name2index[[anctrait]]
    ytree$edge.length <- rep(j, nrow(ytree$edge))     # default every branch = ancestral
    is_term <- ytree$edge[, 2] <= Ntip
    tip_of_edge <- rep(NA_character_, nrow(ytree$edge))
    tip_of_edge[is_term] <- ytree$tip.label[ytree$edge[is_term, 2]]
    for (cat in categories) {
      if (cat == anctrait) next
      sel <- !is.na(tip_of_edge) & pheno[tip_of_edge] == cat
      ytree$edge.length[sel] <- name2index[[cat]]
    }
  }

  ytree$edge.length <- as.numeric(ytree$edge.length)
  list(tree = ytree, legend = legend, n_states = k)
}

# ---------------------------------------------------------------------------
# 3. Write the categorical phenotype tree in CAFE-5 "-y" (lambda) Newick format
# ---------------------------------------------------------------------------
#
# Emits integers (the phenotype class codes) in the branch-length slot, with no
# internal node labels and no root edge, e.g.
#     ((cat:2,horse:2):2,(rat:1,mouse:1):1);
# `tree` is normally the $tree returned by build_categorical_phenotype_tree().
# Returns `path` invisibly.
write_categorical_phenotype_tree <- function(tree, path) {
  .cpt_require_pkg("ape")
  if (is.null(tree$edge.length)) {
    stop("`tree` has no edge.length (the integer class codes).", call. = FALSE)
  }
  out <- tree
  out$node.label <- NULL
  out$root.edge  <- NULL
  out$edge.length <- as.integer(round(out$edge.length))
  # write.tree writes integer edge lengths cleanly; strip any accidental
  # trailing decimals defensively so CAFE sees ":2" not ":2.0000000".
  nwk <- ape::write.tree(out, file = "")
  nwk <- gsub(":([0-9]+)\\.0+([,)])", ":\\1\\2", nwk)
  writeLines(nwk, path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# 4. Prune a categorical phenotype tree without corrupting the integer codes
# ---------------------------------------------------------------------------
#
# Uses categoricalDropTip() (from CategoricalDropTip.R), which -- unlike
# ape::drop.tip -- keeps the CHILD's value when collapsing single-child nodes
# instead of SUMMING branch lengths (see CategoricalDropTip.R). Summing would be
# meaningless for integer category codes.
#
# Use ONLY when the phenotyped set is a strict subset of the tree's tips AND you
# also prune the "-t" tree and the CAFE count matrix to the same tip set, so all
# three CAFE inputs stay consistent.
prune_categorical_phenotype_tree <- function(tree, drop_tips) {
  if (!exists("categoricalDropTip", mode = "function")) {
    stop(paste0("categoricalDropTip() not found; source ",
                "Src/Reu/r/CategoricalDropTip.R before calling ",
                "prune_categorical_phenotype_tree()."), call. = FALSE)
  }
  categoricalDropTip(tree, drop_tips)
}
