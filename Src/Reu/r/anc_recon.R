# anc_recon.R
#
# VENDORED ancestral-state-reconstruction engine for the categorical phenotype
# tree. This is a copy of getAncLiks() from the Meyer lab's RERconverge code:
#     wkmeyer-lab/RunRERBinaryMT @ Dev : Src/Reu/RERConvergeFunctions.R
# reproduced here so the categorical-phenotype-tree step depends on a couple of
# small CRAN packages instead of the whole RERconverge package. RERconverge's
# char2TreeCategorical() calls this exact routine internally, so the states it
# produces match the lab's pipeline (the "Tier-2 exact-parity" option).
#
# WHAT getAncLiks() DOES -- marginal ancestral-state reconstruction of a discrete
# trait under a continuous-time Markov (Mk) model:
#   1. one-hot the observed tip states,
#   2. fit the k x k transition-rate matrix Q with castor::fit_mk
#      (rate_model = "ER" | "SYM" | "ARD"),
#   3. Felsenstein pruning FORWARD pass (postorder): propagate conditional
#      likelihoods up the tree with P = expm(Q * branch_length),
#   4. BACKWARD pass (preorder): combine each node with the rest of the tree to
#      get its marginal posterior over the k states.
# Returns a matrix with one row per INTERNAL node (ordered by node index
# Ntip+1, Ntip+2, ...) and one column per state (states 1..k in sorted order).
#
# DEPENDENCIES (small; installable from CRAN; NO RERconverge):
#   - castor : fit_mk()  -- the Mk-model fit; the numerical heart.
#   - expm   : expm()    -- matrix exponential; kept as-is for exact parity.
#   - ape (namespace)    -- reorder() dispatch for phylo objects.
#   (Matrix / phangorn are NOT needed; see the substitutions below.)
#
# DELIBERATE, PARITY-PRESERVING SUBSTITUTIONS vs. the lab source:
#   - the lab builds the tip one-hot with `to.matrix(tipvals, sort(unique(...)))`;
#     we build the identical 0/1 matrix in base R with outer(), so phangorn is
#     not pulled in. The matrix is bit-for-bit the same.
#   - the lab transposes a dense vector with `Matrix::t`; we use base t(), which
#     is identical for the dense base matrices produced here, so Matrix is not
#     required as an explicit dependency.
#   - castor::fit_mk and expm::expm are called with explicit `::` so the packages
#     only need to be installed (requireNamespace), not attached.
#
# tipvals MUST be integer state codes (1..k) ordered to tree$tip.label (nodes
# 1..Ntip); this is how build_categorical_phenotype_tree() calls it.

getAncLiks <- function(tree, tipvals, Q = NULL, rate_model = "ER", root_prior = "auto") {
  for (pkg in c("ape", "castor", "expm")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required by the vendored getAncLiks().", pkg),
           call. = FALSE)
    }
  }

  ntips <- length(tree$tip.label)
  states_sorted <- sort(unique(tipvals))
  nstates <- length(states_sorted)

  # --- tip one-hot likelihoods (base-R replacement for phangorn::to.matrix) ---
  # tips[i, j] = 1 if tip i is in state states_sorted[j], else 0.
  tips <- outer(tipvals, states_sorted, FUN = "==") * 1
  liks <- matrix(nrow = tree$Nnode, ncol = nstates)
  liks <- rbind(tips, liks)   # rows 1..ntips = tips; ntips+1.. = internal nodes

  # --- fit the Mk transition-rate matrix Q ---
  if (is.null(Q)) {
    Q <- castor::fit_mk(
      trees = tree, Nstates = nstates, tip_states = tipvals,
      rate_model = rate_model, root_prior = root_prior
    )$transition_matrix
  }

  tree <- reorder(tree, order = "postorder")   # ape's reorder.phylo (namespace-loaded)

  # --- forward pass (Felsenstein pruning) ---
  parents <- unique(tree$edge[, 1])
  for (i in 1:length(parents)) {
    p <- parents[i]
    cc <- tree$edge[, 2][which(tree$edge[, 1] == p)]
    ee <- tree$edge.length[which(tree$edge[, 1] == p)]
    v <- vector(mode = "list", length = length(cc))
    for (c in 1:length(cc)) {
      P <- expm::expm(Q * ee[c])
      v[[c]] <- P %*% liks[cc[c], ]
    }
    ll <- Reduce("*", v)[, 1]
    liks[p, ] <- ll / sum(ll)
  }

  # --- backward pass (marginal reconstruction) ---
  for (i in length(parents):1) {
    p <- parents[i]
    cc <- tree$edge[, 2][which(tree$edge[, 1] == p)]
    ee <- tree$edge.length[which(tree$edge[, 1] == p)]
    for (c in 1:length(cc)) {
      des <- cc[c]
      if (des > ntips) {
        P <- expm::expm(Q * ee[c])
        tmp <- t(liks[p, ] / (P %*% liks[des, ]))
        ll <- (tmp %*% P) * liks[des, ]
        ll <- ll[1, ]
        liks[des, ] <- ll / sum(ll, na.rm = TRUE)
      }
    }
  }

  liks <- liks[-(1:ntips), ]
  if (sum(is.nan(liks)) > 0) {
    warning(paste0("NaN values produced. Zero values in the transition rate ",
                   "matrix preventing necessary transitions. Use a different ",
                   "rate model."), call. = FALSE)
  }
  liks
}
