#' Phylogenetic covariance models
#'
#' Functions to compute phylogenetic covariance matrices under Brownian motion
#' (BM), Pagel's lambda-transformed BM, and Ornstein-Uhlenbeck (OU) models.
#'
#' @name covariance
NULL

#' @rdname covariance
#' @export
#' @title Compute phylogenetic covariance matrix
#' @description Computes a phylogenetic covariance matrix under one of three
#'   models: Brownian motion (BM), Pagel's lambda-transformed BM, or
#'   Ornstein-Uhlenbeck (OU).
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param tips Optional character vector of tip names to include. If NULL,
#'   all tips are used.
#' @param model Covariance model: \code{"BM"}, \code{"lambda"}, or \code{"OU"}.
#' @param lambda Lambda parameter for Pagel's lambda model (0 to 1).
#'   Ignored unless \code{model = "lambda"}.
#' @param half_life Half-life parameter for OU model. If provided, overrides
#'   \code{alpha}. Ignored unless \code{model = "OU"}.
#' @param alpha OU selection strength parameter. If \code{half_life} is
#'   provided, alpha is computed as log(2) / half_life.
#' @return A numeric covariance matrix with tip names as row and column names.
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(10)
#' V <- phylo_covariance(tree, model = "BM")
#' dim(V)
#' }
phylo_covariance <- function(tree, tips = NULL,
                             model = c("BM", "lambda", "OU"),
                             lambda = 1, half_life = NULL, alpha = NULL) {
  model <- match.arg(model)

  if (model == "BM") {
    V <- ape::vcv.phylo(tree, corr = FALSE)
  } else if (model == "lambda") {
    V_bm <- ape::vcv.phylo(tree, corr = FALSE)
    V <- lambda_transform_cov(V_bm, lambda)
  } else if (model == "OU") {
    if (!is.null(half_life)) {
      alpha <- log(2) / half_life
    }
    if (is.null(alpha) || alpha <= 0) {
      stop("For OU model, provide either alpha (> 0) or half_life (> 0)")
    }
    V <- make_ou_covariance(tree, alpha)
  }

  if (!is.null(tips)) {
    missing <- setdiff(tips, rownames(V))
    if (length(missing) > 0) {
      stop(paste0(length(missing), " specified tips not found in covariance matrix"))
    }
    V <- V[tips, tips, drop = FALSE]
  }

  V
}

#' Lambda-transform a BM covariance matrix
#'
#' Off-diagonal elements are multiplied by lambda; diagonal stays unchanged.
#'
#' @param V_bm BM covariance matrix
#' @param lambda Scalar between 0 and 1
#' @return Transformed covariance matrix
#' @keywords internal
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
#' @param tree A phylo object
#' @param alpha OU selection strength parameter (positive)
#' @return OU covariance matrix
#' @keywords internal
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

#' @rdname covariance
#' @export
#' @title Convert covariance matrix to correlation matrix
#' @description Standardizes a covariance matrix to a correlation matrix.
#' @param V A covariance matrix.
#' @return A correlation matrix.
#' @examples
#' V <- matrix(c(4, 1, 1, 9), nrow = 2)
#' cov_to_cor(V)
cov_to_cor <- function(V) {
  stats::cov2cor(V)
}

#' @rdname covariance
#' @export
#' @title Compute dependence metrics from a correlation matrix
#' @description Computes MeanOffCor (mean absolute off-diagonal correlation),
#'   MaxOffCor (maximum absolute off-diagonal correlation), and MeanESS
#'   (mean-based effective sample size) from a correlation matrix.
#' @param R A correlation matrix.
#' @return A list with components:
#'   \item{MeanOffCor}{Mean absolute off-diagonal correlation.}
#'   \item{MaxOffCor}{Maximum absolute off-diagonal correlation.}
#'   \item{MeanESS}{Mean-based effective sample size.}
#' @examples
#' R <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' dependence_metrics(R)
dependence_metrics <- function(R) {
  R <- as.matrix(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1

  off_values <- abs(R[upper.tri(R)])

  n <- nrow(R)
  one_vec <- rep(1, n)

  ess <- tryCatch({
    chol_R <- chol(R)
    inv_R <- chol2inv(chol_R)
    as.numeric(t(one_vec) %*% inv_R %*% one_vec)
  }, error = function(e) {
    tryCatch({
      as.numeric(t(one_vec) %*% solve(R, one_vec))
    }, error = function(e2) {
      NA_real_
    })
  })

  list(
    MeanOffCor = mean(off_values, na.rm = TRUE),
    MaxOffCor = max(off_values, na.rm = TRUE),
    MeanESS = ess
  )
}

#' @rdname covariance
#' @export
#' @title Mean off-diagonal correlation
#' @param R A correlation matrix.
#' @return Mean absolute off-diagonal correlation.
mean_off_cor <- function(R) {
  R <- as.matrix(R)
  off_values <- abs(R[upper.tri(R)])
  mean(off_values, na.rm = TRUE)
}

#' @rdname covariance
#' @export
#' @title Maximum off-diagonal correlation
#' @param R A correlation matrix.
#' @return Maximum absolute off-diagonal correlation.
max_off_cor <- function(R) {
  R <- as.matrix(R)
  off_values <- abs(R[upper.tri(R)])
  max(off_values, na.rm = TRUE)
}

#' @rdname covariance
#' @export
#' @title Mean-based effective sample size
#' @description Computes the mean-based effective sample size (MeanESS) as
#'   1' R^{-1} 1, where R is the correlation matrix of the subset. This
#'   quantifies the effective information content for estimating a mean or
#'   similar aggregate quantity.
#' @param R A correlation matrix.
#' @param method Method for matrix inversion: \code{"qr"} (default, using
#'   \code{qr.solve}), \code{"chol"} (Cholesky decomposition), or
#'   \code{"ginv"} (pseudo-inverse via \code{ginv}).
#' @param tol Tolerance for detecting singularity.
#' @return The effective sample size (numeric). Returns NA if the matrix is
#'   singular and cannot be inverted.
#' @examples
#' R <- diag(5)
#' mean_ess(R)  # Should be 5
mean_ess <- function(R, method = c("qr", "chol", "ginv"), tol = 1e-8) {
  method <- match.arg(method)
  R <- as.matrix(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1

  n <- nrow(R)
  one_vec <- rep(1, n)

  result <- tryCatch({
    if (method == "qr") {
      inv_R <- qr.solve(R, tol = tol)
      as.numeric(t(one_vec) %*% inv_R %*% one_vec)
    } else if (method == "chol") {
      chol_R <- chol(R)
      inv_R <- chol2inv(chol_R)
      as.numeric(t(one_vec) %*% inv_R %*% one_vec)
    } else {
      # ginv
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("MASS package required for method = 'ginv'")
      }
      inv_R <- MASS::ginv(R, tol = tol)
      as.numeric(t(one_vec) %*% inv_R %*% one_vec)
    }
  }, error = function(e) {
    NA_real_
  })

  result
}
