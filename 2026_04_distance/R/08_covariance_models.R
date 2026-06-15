# 08_covariance_models.R
# Covariance models: BM, lambda-transformed BM, OU, EB
# Borrows OU logic from old signal_metrics.R but implements lambda as
# off-diagonal rescaling (not branch-length transformation).

#' Make BM covariance matrix (vcv)
#'
#' @param tree A phylo object
#' @return BM covariance matrix
make_bm_covariance <- function(tree) {
  ape::vcv.phylo(tree, corr = FALSE)
}

#' Lambda-transform a BM covariance matrix
#'
#' Off-diagonal elements are multiplied by lambda; diagonal stays unchanged.
#' This is the sensitivity scheme described in the manuscript.
#'
#' @param V_bm BM covariance matrix
#' @param lambda Scalar between 0 and 1
#' @return Transformed covariance matrix
lambda_transform_cov <- function(V_bm, lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("lambda must be a single numeric value.")
  }
  if (lambda < 0 || lambda > 1) {
    stop("lambda should be between 0 and 1.")
  }
  
  V_lam <- V_bm
  offdiag <- row(V_lam) != col(V_lam)
  V_lam[offdiag] <- lambda * V_lam[offdiag]
  
  V_lam
}

#' Make OU covariance matrix
#'
#' Based on root-to-tip times, MRCA shared times, and pairwise patristic distances.
#'
#' @param tree A phylo object
#' @param alpha OU selection strength parameter (positive)
#' @return OU covariance matrix
make_ou_covariance <- function(tree, alpha) {
  n_tips <- ape::Ntip(tree)
  
  if (alpha <= 0) {
    stop("alpha must be positive.")
  }
  
  node_depths <- ape::node.depth.edgelength(tree)
  tip_depths <- node_depths[seq_len(n_tips)]
  names(tip_depths) <- tree$tip.label
  
  mrca_matrix <- ape::mrca(tree)
  mrca_tips <- mrca_matrix[tree$tip.label, tree$tip.label, drop = FALSE]
  
  shared_times <- matrix(
    node_depths[mrca_tips],
    nrow = n_tips,
    ncol = n_tips,
    dimnames = list(tree$tip.label, tree$tip.label)
  )
  
  dist_matrix <- ape::cophenetic.phylo(tree)
  
  cov_matrix <- (1 / (2 * alpha)) *
    exp(-alpha * dist_matrix) *
    (1 - exp(-2 * alpha * shared_times))
  
  diag(cov_matrix) <- (1 / (2 * alpha)) *
    (1 - exp(-2 * alpha * tip_depths))
  
  cov_matrix
}

#' Convert half-life to OU alpha
#'
#' @param half_life Half-life value
#' @return alpha = log(2) / half_life
half_life_to_alpha <- function(half_life) {
  log(2) / half_life
}

#' Make OU covariance by half-life fraction of tree height
#'
#' @param tree A phylo object
#' @param half_life_frac Fraction of tree height (e.g. 0.1, 0.25, 0.5, 1.0)
#' @return List with V, alpha, half_life, half_life_frac
make_ou_covariance_by_half_life_fraction <- function(tree, half_life_frac) {
  tree_height <- max(ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))])
  half_life <- half_life_frac * tree_height
  alpha <- half_life_to_alpha(half_life)
  
  list(
    V = make_ou_covariance(tree, alpha),
    alpha = alpha,
    half_life = half_life,
    half_life_frac = half_life_frac
  )
}

#' Make Early-Burst covariance matrix
#'
#' Branch lengths are transformed by integrating exp(rate * t) from the
#' parent node time to the child node time, then passed to ape::vcv.phylo.
#'
#' @param tree A phylo object
#' @param rate Numeric scalar EB rate
#' @return EB covariance matrix
make_eb_covariance <- function(tree, rate) {
  if (!is.numeric(rate) || length(rate) != 1 || !is.finite(rate)) {
    stop("rate must be a single finite numeric value.")
  }
  if (abs(rate) < 1e-8) {
    return(make_bm_covariance(tree))
  }
  if (is.null(tree$edge.length)) {
    stop("tree must have edge lengths.")
  }
  
  node_depths <- ape::node.depth.edgelength(tree)
  parent_times <- node_depths[tree$edge[, 1]]
  child_times <- node_depths[tree$edge[, 2]]
  
  eb_tree <- tree
  eb_tree$edge.length <- (exp(rate * child_times) - exp(rate * parent_times)) / rate
  eb_tree$edge.length <- pmax(eb_tree$edge.length, 0)
  
  V_eb <- ape::vcv.phylo(eb_tree, corr = FALSE)
  V_eb <- V_eb[tree$tip.label, tree$tip.label, drop = FALSE]
  rownames(V_eb) <- tree$tip.label
  colnames(V_eb) <- tree$tip.label
  V_eb
}

#' Make simple EB-style covariance sensitivity matrix
#'
#' Starting from BM covariance, rescale shared-history covariance by
#' exp(-r * t_scaled), where t_scaled is shared ancestry / tree height.
#'
#' @param tree A phylo object
#' @param r Numeric scalar EB sensitivity rate
#' @return EB-style covariance matrix
make_eb_covariance_simple <- function(tree, r) {
  if (!is.numeric(r) || length(r) != 1 || !is.finite(r)) {
    stop("r must be a single finite numeric value.")
  }
  
  V_bm <- make_bm_covariance(tree)
  node_depths <- ape::node.depth.edgelength(tree)
  n_tips <- ape::Ntip(tree)
  tree_height <- max(node_depths[seq_len(n_tips)])
  
  if (!is.finite(tree_height) || tree_height <= 0) {
    stop("tree height must be positive and finite.")
  }
  
  mrca_mat <- ape::mrca(tree)
  mrca_tips <- mrca_mat[tree$tip.label, tree$tip.label, drop = FALSE]
  
  shared_time <- matrix(
    node_depths[mrca_tips],
    nrow = n_tips,
    ncol = n_tips,
    dimnames = list(tree$tip.label, tree$tip.label)
  )
  
  t_scaled <- shared_time / tree_height
  V_eb <- V_bm * exp(-r * t_scaled)
  V_eb <- (V_eb + t(V_eb)) / 2
  rownames(V_eb) <- tree$tip.label
  colnames(V_eb) <- tree$tip.label
  V_eb
}

