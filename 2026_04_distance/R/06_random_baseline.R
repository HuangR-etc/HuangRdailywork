# 06_random_baseline.R
# Random subset sampling, p-value and SES calculation
# Migrated from old subset_random.R

#' Sample random subsets from a tree
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @param n_reps Number of random subsets to sample
#' @param replace If TRUE, sample with replacement (default: FALSE)
#' @return A list of random subsets (each as vector of indices)
sample_random_subsets <- function(dist_obj, subset_size, n_reps, replace = FALSE) {
  all_tips <- seq_along(dist_obj$tip_labels)
  
  if (subset_size > length(all_tips)) {
    stop("Subset size cannot be larger than total number of tips")
  }
  
  random_subsets <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    if (replace) {
      random_subsets[[i]] <- sample(all_tips, subset_size, replace = TRUE)
    } else {
      random_subsets[[i]] <- sample(all_tips, subset_size, replace = FALSE)
    }
  }
  
  random_subsets
}

#' Calculate one-sided p-value for observed >= null (upper tail)
#'
#' @param obs Observed value
#' @param null_values Numeric vector of null distribution values
#' @return P-value
calc_p_high <- function(obs, null_values) {
  (1 + sum(null_values >= obs)) / (1 + length(null_values))
}

#' Calculate one-sided p-value for observed <= null (lower tail)
#'
#' @param obs Observed value
#' @param null_values Numeric vector of null distribution values
#' @return P-value
calc_p_low <- function(obs, null_values) {
  (1 + sum(null_values <= obs)) / (1 + length(null_values))
}

#' Calculate Standardized Effect Size (SES)
#'
#' @param obs Observed value
#' @param null_values Numeric vector of null distribution values
#' @return SES (z-score), NA if sd = 0
calc_ses <- function(obs, null_values) {
  sd_x <- sd(null_values)
  if (!is.finite(sd_x) || sd_x == 0) return(NA_real_)
  (obs - mean(null_values)) / sd_x
}

#' Calculate null distribution metrics for random subsets (extended)
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @param n_reps Number of random subsets
#' @return A data frame with extended metrics for all random subsets
calc_null_distribution_extended <- function(dist_obj, subset_size, n_reps) {
  random_subsets <- sample_random_subsets(dist_obj, subset_size, n_reps)
  calc_multiple_subsets_metrics_extended(dist_obj$dist_mat, random_subsets)
}
