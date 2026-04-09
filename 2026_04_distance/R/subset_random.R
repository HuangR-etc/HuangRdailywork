# Random subset sampling and null distribution analysis
# This module handles random subset sampling and empirical p-value calculation

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
  
  random_subsets <- list()
  
  for (i in 1:n_reps) {
    if (replace) {
      # Sample with replacement (not typical for subset selection)
      random_subsets[[i]] <- sample(all_tips, subset_size, replace = TRUE)
    } else {
      # Sample without replacement (typical)
      random_subsets[[i]] <- sample(all_tips, subset_size, replace = FALSE)
    }
  }
  
  return(random_subsets)
}

#' Calculate null distribution metrics for random subsets
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @param n_reps Number of random subsets
#' @return A data frame with metrics for all random subsets
calc_null_distribution <- function(dist_obj, subset_size, n_reps) {
  cat("Calculating null distribution with", n_reps, "random subsets...\n")
  
  # Sample random subsets
  random_subsets <- sample_random_subsets(dist_obj, subset_size, n_reps)
  
  # Calculate metrics for each subset
  null_metrics <- calc_multiple_subsets_metrics(dist_obj$dist_mat, random_subsets)
  
  # Add subset information
  null_metrics$SubsetType <- "random"
  null_metrics$SubsetSize <- subset_size
  
  cat("Null distribution calculation complete.\n")
  return(null_metrics)
}

#' Calculate empirical p-value for observed metrics
#'
#' @param observed_metrics List with MinPD, MeanPD, MeanNND for observed subset
#' @param null_metrics Data frame with null distribution metrics
#' @param maximize If TRUE, test if observed is larger than null (default: TRUE)
#' @return A list with p-values for each metric
empirical_pvalue <- function(observed_metrics, null_metrics, maximize = TRUE) {
  n_null <- nrow(null_metrics)
  
  if (maximize) {
    # For maximization: p = (1 + #(null >= observed)) / (1 + n_null)
    p_minpd <- (1 + sum(null_metrics$MinPD >= observed_metrics$MinPD)) / (1 + n_null)
    p_meanpd <- (1 + sum(null_metrics$MeanPD >= observed_metrics$MeanPD)) / (1 + n_null)
    p_meannnd <- (1 + sum(null_metrics$MeanNND >= observed_metrics$MeanNND)) / (1 + n_null)
  } else {
    # For minimization: p = (1 + #(null <= observed)) / (1 + n_null)
    p_minpd <- (1 + sum(null_metrics$MinPD <= observed_metrics$MinPD)) / (1 + n_null)
    p_meanpd <- (1 + sum(null_metrics$MeanPD <= observed_metrics$MeanPD)) / (1 + n_null)
    p_meannnd <- (1 + sum(null_metrics$MeanNND <= observed_metrics$MeanNND)) / (1 + n_null)
  }
  
  return(list(
    MinPD = p_minpd,
    MeanPD = p_meanpd,
    MeanNND = p_meannnd
  ))
}

#' Calculate z-score for observed metrics relative to null distribution
#'
#' @param observed_metrics List with MinPD, MeanPD, MeanNND for observed subset
#' @param null_metrics Data frame with null distribution metrics
#' @return A list with z-scores for each metric
calculate_zscore <- function(observed_metrics, null_metrics) {
  z_minpd <- (observed_metrics$MinPD - mean(null_metrics$MinPD)) / sd(null_metrics$MinPD)
  z_meanpd <- (observed_metrics$MeanPD - mean(null_metrics$MeanPD)) / sd(null_metrics$MeanPD)
  z_meannnd <- (observed_metrics$MeanNND - mean(null_metrics$MeanNND)) / sd(null_metrics$MeanNND)
  
  return(list(
    MinPD = z_minpd,
    MeanPD = z_meanpd,
    MeanNND = z_meannnd
  ))
}

#' Calculate percentile rank for observed metrics
#'
#' @param observed_metrics List with MinPD, MeanPD, MeanNND for observed subset
#' @param null_metrics Data frame with null distribution metrics
#' @param maximize If TRUE, calculate percentile for maximization (default: TRUE)
#' @return A list with percentiles for each metric
calculate_percentile <- function(observed_metrics, null_metrics, maximize = TRUE) {
  if (maximize) {
    # Percentile = proportion of null values <= observed
    pct_minpd <- mean(null_metrics$MinPD <= observed_metrics$MinPD)
    pct_meanpd <- mean(null_metrics$MeanPD <= observed_metrics$MeanPD)
    pct_meannnd <- mean(null_metrics$MeanNND <= observed_metrics$MeanNND)
  } else {
    # For minimization: percentile = proportion of null values >= observed
    pct_minpd <- mean(null_metrics$MinPD >= observed_metrics$MinPD)
    pct_meanpd <- mean(null_metrics$MeanPD >= observed_metrics$MeanPD)
    pct_meannnd <- mean(null_metrics$MeanNND >= observed_metrics$MeanNND)
  }
  
  return(list(
    MinPD = pct_minpd,
    MeanPD = pct_meanpd,
    MeanNND = pct_meannnd
  ))
}

#' Compare observed subset with null distribution
#'
#' @param dist_obj Distance object
#' @param observed_subset Observed subset (indices)
#' @param subset_size Subset size
#' @param n_reps Number of random subsets for null distribution
#' @param maximize If TRUE, test if observed is better than random (default: TRUE)
#' @return A comprehensive comparison result
compare_with_null <- function(dist_obj, observed_subset, subset_size, n_reps, maximize = TRUE) {
  # Calculate observed metrics
  observed_metrics <- calc_subset_metrics(dist_obj$dist_mat, observed_subset)
  
  # Calculate null distribution
  null_metrics <- calc_null_distribution(dist_obj, subset_size, n_reps)
  
  # Calculate statistics
  p_values <- empirical_pvalue(observed_metrics, null_metrics, maximize)
  z_scores <- calculate_zscore(observed_metrics, null_metrics)
  percentiles <- calculate_percentile(observed_metrics, null_metrics, maximize)
  
  # Create summary
  summary_df <- data.frame(
    Metric = c("MinPD", "MeanPD", "MeanNND"),
    Observed = c(observed_metrics$MinPD, observed_metrics$MeanPD, observed_metrics$MeanNND),
    Null_Mean = c(mean(null_metrics$MinPD), mean(null_metrics$MeanPD), mean(null_metrics$MeanNND)),
    Null_SD = c(sd(null_metrics$MinPD), sd(null_metrics$MeanPD), sd(null_metrics$MeanNND)),
    Z_Score = c(z_scores$MinPD, z_scores$MeanPD, z_scores$MeanNND),
    P_Value = c(p_values$MinPD, p_values$MeanPD, p_values$MeanNND),
    Percentile = c(percentiles$MinPD, percentiles$MeanPD, percentiles$MeanNND),
    stringsAsFactors = FALSE
  )
  
  return(list(
    observed_metrics = observed_metrics,
    null_metrics = null_metrics,
    summary = summary_df,
    p_values = p_values,
    z_scores = z_scores,
    percentiles = percentiles,
    comparison_type = ifelse(maximize, "maximization", "minimization")
  ))
}

#' Test random subset functions
test_random <- function() {
  # Load required libraries and functions
  library(ape)
  source("distance_metrics.R")
  source("subset_greedy.R")
  
  # Create a test tree
  test_tree <- rtree(30)
  test_tree$tip.label <- paste0("sp", 1:30)
  
  # Create distance object
  dist_obj <- create_distance_object(test_tree)
  
  cat("Testing random subset functions:\n")
  
  # Test 1: Sample random subsets
  cat("\n1. Sampling random subsets:\n")
  random_subsets <- sample_random_subsets(dist_obj, subset_size = 5, n_reps = 10)
  cat("   Sampled", length(random_subsets), "random subsets of size 5\n")
  cat("   First subset indices:", random_subsets[[1]], "\n")
  cat("   First subset names:", dist_obj$tip_labels[random_subsets[[1]]], "\n")
  
  # Test 2: Calculate null distribution
  cat("\n2. Calculating null distribution:\n")
  null_metrics <- calc_null_distribution(dist_obj, subset_size = 5, n_reps = 100)
  cat("   Null distribution statistics:\n")
  cat("     MinPD: mean =", mean(null_metrics$MinPD), "sd =", sd(null_metrics$MinPD), "\n")
  cat("     MeanPD: mean =", mean(null_metrics$MeanPD), "sd =", sd(null_metrics$MeanPD), "\n")
  cat("     MeanNND: mean =", mean(null_metrics$MeanNND), "sd =", sd(null_metrics$MeanNND), "\n")
  
  # Test 3: Create an observed subset (using greedy)
  cat("\n3. Creating observed subset (using greedy):\n")
  greedy_result <- build_subset_greedy(dist_obj, subset_size = 5, maximize = TRUE)
  observed_metrics <- greedy_result$metrics
  cat("   Observed metrics: MinPD =", observed_metrics$MinPD,
      "MeanPD =", observed_metrics$MeanPD,
      "MeanNND =", observed_metrics$MeanNND, "\n")
  
  # Test 4: Compare with null
  cat("\n4. Comparing observed with null distribution:\n")
  comparison <- compare_with_null(dist_obj, greedy_result$subset, 
                                  subset_size = 5, n_reps = 100, maximize = TRUE)
  
  cat("   Comparison summary:\n")
  print(comparison$summary)
  
  # Test 5: Empirical p-values
  cat("\n5. Empirical p-values:\n")
  p_vals <- empirical_pvalue(observed_metrics, null_metrics, maximize = TRUE)
  cat("   MinPD p-value:", p_vals$MinPD, "\n")
  cat("   MeanPD p-value:", p_vals$MeanPD, "\n")
  cat("   MeanNND p-value:", p_vals$MeanNND, "\n")
  
  return(list(
    random_subsets = random_subsets,
    null_metrics = null_metrics,
    comparison = comparison,
    p_values = p_vals
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_random()
}
