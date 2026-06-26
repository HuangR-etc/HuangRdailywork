#' Covariance and dependence diagnostics
#'
#' Implements the BM, lambda-transformed BM, OU, and EB covariance models used
#' in the manuscript, together with MeanOffCor, MaxOffCor, and MIESS.
NULL

CLUSTERED_MAX_EXCHANGE_ITERATIONS <- 10
CLUSTERED_EXCHANGE_TOL <- 1e-10

make_bm_covariance <- function(tree) {
  ape::vcv.phylo(tree, corr = FALSE)
}

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

half_life_to_alpha <- function(half_life) {
  log(2) / half_life
}

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

# Main manuscript analyses used the EB-style shared-history rescaling below.
# Keep the branch-transformation implementation available for backward
# compatibility and method comparisons, but make the manuscript-matching
# version the package default through phylo_covariance().
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

resolve_eb_covariance <- function(tree, rate,
                                  eb_method = c("simple", "branch_transform")) {
  eb_method <- match.arg(eb_method)

  if (eb_method == "simple") {
    return(make_eb_covariance_simple(tree, r = rate))
  }

  make_eb_covariance(tree, rate = rate)
}

#' @export
phylo_covariance <- function(tree, tips = NULL,
                             model = c("BM", "lambda", "OU", "EB"),
                             lambda = 1, half_life = NULL, alpha = NULL,
                             eb_rate = NULL,
                             eb_method = c("simple", "branch_transform")) {
  check_phylo_input(tree)
  model <- match.arg(model)
  eb_method <- match.arg(eb_method)

  V <- switch(
    model,
    BM = make_bm_covariance(tree),
    lambda = lambda_transform_cov(make_bm_covariance(tree), lambda),
    OU = {
      if (!is.null(half_life)) alpha <- log(2) / half_life
      make_ou_covariance(tree, alpha = alpha)
    },
    EB = {
      if (is.null(eb_rate)) {
        stop("eb_rate must be supplied for model = 'EB'.", call. = FALSE)
      }
      resolve_eb_covariance(tree, rate = eb_rate, eb_method = eb_method)
    }
  )

  if (!is.null(tips)) {
    missing_tips <- setdiff(tips, rownames(V))
    if (length(missing_tips) > 0) {
      stop(
        paste0(length(missing_tips), " specified tips not found in covariance matrix."),
        call. = FALSE
      )
    }
    V <- V[tips, tips, drop = FALSE]
  }

  V
}

#' @export
cov_to_cor <- function(V) stats::cov2cor(V)

#' @export
mean_off_cor <- function(R) {
  R <- as.matrix(R)
  mean(R[upper.tri(R)], na.rm = TRUE)
}

#' @export
max_off_cor <- function(R) {
  R <- as.matrix(R)
  max(R[upper.tri(R)], na.rm = TRUE)
}

#' @export
miess <- function(R) {
  R <- as.matrix(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  one <- rep(1, nrow(R))

  tryCatch({
    chol_R <- chol(R)
    inv_R <- chol2inv(chol_R)
    as.numeric(t(one) %*% inv_R %*% one)
  }, error = function(e) {
    tryCatch(as.numeric(t(one) %*% solve(R, one)), error = function(e2) NA_real_)
  })
}

#' @export
dependence_metrics <- function(R) {
  list(
    MeanOffCor = mean_off_cor(R),
    MaxOffCor = max_off_cor(R),
    MIESS = miess(R)
  )
}
