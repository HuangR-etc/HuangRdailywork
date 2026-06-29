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
#' Implements the S4.4 EB construction on a fixed ultrametric BM baseline tree.
#' Branch-specific rates are r(u) = exp(rate * u) on normalized tree time
#' u = t / H, where H is total tree height. Each branch length is replaced by
#' the integral of this rate over original time before passing the transformed
#' tree to ape::vcv.phylo.
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
  if (!isTRUE(ape::is.rooted(tree))) {
    stop("tree must be rooted for EB covariance construction.")
  }
  
  node_depths <- ape::node.depth.edgelength(tree)
  tip_depths <- node_depths[seq_len(ape::Ntip(tree))]
  ultrametric_tol <- 1e-6 * max(1, max(tip_depths))
  if ((max(tip_depths) - min(tip_depths)) > ultrametric_tol) {
    stop("tree must be ultrametric for EB covariance construction.")
  }
  tree_height <- mean(tip_depths)
  if (!is.finite(tree_height) || tree_height <= 0) {
    stop("tree height must be positive for EB covariance construction.")
  }
  parent_times <- node_depths[tree$edge[, 1]] / tree_height
  child_times <- node_depths[tree$edge[, 2]] / tree_height
  
  eb_tree <- tree
  eb_tree$edge.length <- tree_height *
    (exp(rate * child_times) - exp(rate * parent_times)) / rate
  eb_tree$edge.length <- pmax(eb_tree$edge.length, 0)
  
  V_eb <- ape::vcv.phylo(eb_tree, corr = FALSE)
  V_eb <- V_eb[tree$tip.label, tree$tip.label, drop = FALSE]
  rownames(V_eb) <- tree$tip.label
  colnames(V_eb) <- tree$tip.label
  V_eb
}
