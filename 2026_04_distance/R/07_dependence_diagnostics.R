# 07_dependence_diagnostics.R
# Dependence diagnostics from covariance/correlation matrices
# Calculates off_mean, rmax, neff_mean (MeanESS = 1' R^{-1} 1)

#' Calculate dependence diagnostics from a correlation matrix
#'
#' @param R_sub Correlation matrix for the subset
#' @return Data frame with off_mean, rmax, neff_mean
calc_dependence_from_R <- function(R_sub) {
  R_sub <- as.matrix(R_sub)
  R_sub <- (R_sub + t(R_sub)) / 2
  diag(R_sub) <- 1
  
  off_values <- R_sub[upper.tri(R_sub)]
  
  n <- nrow(R_sub)
  one_vec <- rep(1, n)
  
  neff_mean <- tryCatch({
    chol_R <- chol(R_sub)
    inv_R <- chol2inv(chol_R)
    as.numeric(t(one_vec) %*% inv_R %*% one_vec)
  }, error = function(e) {
    tryCatch({
      as.numeric(t(one_vec) %*% solve(R_sub, one_vec))
    }, error = function(e2) {
      NA_real_
    })
  })
  
  data.frame(
    off_mean = mean(off_values, na.rm = TRUE),
    rmax = max(off_values, na.rm = TRUE),
    neff_mean = neff_mean
  )
}

#' Calculate dependence diagnostics from a full covariance matrix
#'
#' Extracts the subset rows/cols, converts to correlation, then computes diagnostics.
#'
#' @param V_full Full covariance matrix (named)
#' @param subset_names Character vector of tip names in the subset
#' @return Data frame with off_mean, rmax, neff_mean
calc_dependence_from_V <- function(V_full, subset_names) {
  V_sub <- V_full[subset_names, subset_names, drop = FALSE]
  R_sub <- cov2cor(V_sub)
  calc_dependence_from_R(R_sub)
}

#' Calculate dependence diagnostics for multiple subsets
#'
#' @param V_full Full covariance matrix (named)
#' @param subset_names_list List of character vectors (each is tip names for a subset)
#' @return Data frame with off_mean, rmax, neff_mean, SubsetID
calc_multiple_dependence_from_V <- function(V_full, subset_names_list) {
  out <- do.call(rbind, lapply(seq_along(subset_names_list), function(i) {
    x <- subset_names_list[[i]]
    res <- calc_dependence_from_V(V_full, x)
    res$SubsetID <- i
    res
  }))
  out
}
