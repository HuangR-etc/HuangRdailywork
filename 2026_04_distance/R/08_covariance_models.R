# 08_covariance_models.R
# Covariance models: BM, lambda-transformed BM, OU
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
