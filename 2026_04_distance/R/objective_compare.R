# Objective comparison functions for phylogenetic dispersed subset analysis
# This module implements lexicographic comparison of subsets

#' Compare two subsets using lexicographic ordering (maximization)
#'
#' @param metrics_a List with MinPD, MeanPD, MeanNND for subset A
#' @param metrics_b List with MinPD, MeanPD, MeanNND for subset B
#' @return TRUE if subset A is better than subset B, FALSE otherwise
is_better_lexico_max <- function(metrics_a, metrics_b) {
  # Compare MinPD first
  if (metrics_a$MinPD > metrics_b$MinPD) {
    return(TRUE)
  }
  if (metrics_a$MinPD < metrics_b$MinPD) {
    return(FALSE)
  }
  
  # MinPD equal, compare MeanPD
  if (metrics_a$MeanPD > metrics_b$MeanPD) {
    return(TRUE)
  }
  if (metrics_a$MeanPD < metrics_b$MeanPD) {
    return(FALSE)
  }
  
  # MinPD and MeanPD equal, compare MeanNND
  if (metrics_a$MeanNND > metrics_b$MeanNND) {
    return(TRUE)
  }
  
  # Otherwise, B is at least as good as A
  return(FALSE)
}

#' Compare two subsets using lexicographic ordering (minimization)
#'
#' @param metrics_a List with MinPD, MeanPD, MeanNND for subset A
#' @param metrics_b List with MinPD, MeanPD, MeanNND for subset B
#' @return TRUE if subset A is better than subset B, FALSE otherwise
is_better_lexico_min <- function(metrics_a, metrics_b) {
  # Compare MinPD first (minimization)
  if (metrics_a$MinPD < metrics_b$MinPD) {
    return(TRUE)
  }
  if (metrics_a$MinPD > metrics_b$MinPD) {
    return(FALSE)
  }
  
  # MinPD equal, compare MeanPD (minimization)
  if (metrics_a$MeanPD < metrics_b$MeanPD) {
    return(TRUE)
  }
  if (metrics_a$MeanPD > metrics_b$MeanPD) {
    return(FALSE)
  }
  
  # MinPD and MeanPD equal, compare MeanNND (minimization)
  if (metrics_a$MeanNND < metrics_b$MeanNND) {
    return(TRUE)
  }
  
  # Otherwise, B is at least as good as A
  return(FALSE)
}

#' Rank subsets using lexicographic ordering
#'
#' @param subsets List of subsets (each is a vector of tip indices/names)
#' @param dist_mat Distance matrix
#' @param maximize If TRUE, rank by maximization; if FALSE, rank by minimization
#' @return A data frame with ranked subsets and their metrics
rank_subsets_lexico <- function(subsets, dist_mat, maximize = TRUE) {
  # Calculate metrics for all subsets
  metrics_list <- lapply(subsets, function(subset) {
    calc_subset_metrics(dist_mat, subset)
  })
  
  # Create a data frame for sorting
  df <- data.frame(
    SubsetID = seq_along(subsets),
    MinPD = sapply(metrics_list, function(x) x$MinPD),
    MeanPD = sapply(metrics_list, function(x) x$MeanPD),
    MeanNND = sapply(metrics_list, function(x) x$MeanNND),
    stringsAsFactors = FALSE
  )
  
  # Sort based on maximize parameter
  if (maximize) {
    # Sort by MinPD (desc), then MeanPD (desc), then MeanNND (desc)
    df <- df[order(-df$MinPD, -df$MeanPD, -df$MeanNND), ]
  } else {
    # Sort by MinPD (asc), then MeanPD (asc), then MeanNND (asc)
    df <- df[order(df$MinPD, df$MeanPD, df$MeanNND), ]
  }
  
  # Add rank
  df$Rank <- seq_len(nrow(df))
  
  # Reorder columns
  df <- df[, c("Rank", "SubsetID", "MinPD", "MeanPD", "MeanNND")]
  
  return(df)
}

#' Find the best subset from a list using lexicographic ordering
#'
#' @param subsets List of subsets (each is a vector of tip indices/names)
#' @param dist_mat Distance matrix
#' @param maximize If TRUE, find best by maximization; if FALSE, find best by minimization
#' @return The best subset (as vector of indices/names)
find_best_subset <- function(subsets, dist_mat, maximize = TRUE) {
  if (length(subsets) == 0) {
    return(NULL)
  }
  
  # Calculate metrics for all subsets
  metrics_list <- lapply(subsets, function(subset) {
    calc_subset_metrics(dist_mat, subset)
  })
  
  # Initialize with first subset
  best_idx <- 1
  best_metrics <- metrics_list[[1]]
  
  # Compare with remaining subsets
  for (i in 2:length(subsets)) {
    current_metrics <- metrics_list[[i]]
    
    if (maximize) {
      if (is_better_lexico_max(current_metrics, best_metrics)) {
        best_idx <- i
        best_metrics <- current_metrics
      }
    } else {
      if (is_better_lexico_min(current_metrics, best_metrics)) {
        best_idx <- i
        best_metrics <- current_metrics
      }
    }
  }
  
  return(subsets[[best_idx]])
}

#' Compare two subsets directly (wrapper function)
#'
#' @param subset_a First subset (vector of indices/names)
#' @param subset_b Second subset (vector of indices/names)
#' @param dist_mat Distance matrix
#' @param maximize If TRUE, compare for maximization; if FALSE, compare for minimization
#' @return 1 if subset_a is better, -1 if subset_b is better, 0 if equal
compare_subsets <- function(subset_a, subset_b, dist_mat, maximize = TRUE) {
  metrics_a <- calc_subset_metrics(dist_mat, subset_a)
  metrics_b <- calc_subset_metrics(dist_mat, subset_b)
  
  if (maximize) {
    if (is_better_lexico_max(metrics_a, metrics_b)) {
      return(1)
    } else if (is_better_lexico_max(metrics_b, metrics_a)) {
      return(-1)
    } else {
      return(0)
    }
  } else {
    if (is_better_lexico_min(metrics_a, metrics_b)) {
      return(1)
    } else if (is_better_lexico_min(metrics_b, metrics_a)) {
      return(-1)
    } else {
      return(0)
    }
  }
}

# Test function
if (sys.nframe() == 0) {
  # Load required libraries
  library(ape)
  source("distance_metrics.R")
  
  # Create a test tree
  test_tree <- rtree(10)
  test_tree$tip.label <- paste0("sp", 1:10)
  
  # Create distance matrix
  dist_mat <- calc_distance_matrix(test_tree)
  
  # Create test subsets
  subsets <- list(
    c(1, 3, 5, 7),  # Subset 1
    c(2, 4, 6, 8),  # Subset 2
    c(1, 5, 9, 10), # Subset 3
    c(3, 6, 8, 10)  # Subset 4
  )
  
  # Test lexicographic comparison (maximization)
  cat("Testing lexicographic comparison (maximization):\n")
  
  # Compare subset 1 vs subset 2
  metrics1 <- calc_subset_metrics(dist_mat, subsets[[1]])
  metrics2 <- calc_subset_metrics(dist_mat, subsets[[2]])
  
  cat("Subset 1 metrics: MinPD =", metrics1$MinPD, 
      "MeanPD =", metrics1$MeanPD, 
      "MeanNND =", metrics1$MeanNND, "\n")
  cat("Subset 2 metrics: MinPD =", metrics2$MinPD, 
      "MeanPD =", metrics2$MeanPD, 
      "MeanNND =", metrics2$MeanNND, "\n")
  
  result <- is_better_lexico_max(metrics1, metrics2)
  cat("Is subset 1 better than subset 2?", result, "\n")
  
  # Test ranking
  cat("\nRanking subsets (maximization):\n")
  ranked <- rank_subsets_lexico(subsets, dist_mat, maximize = TRUE)
  print(ranked)
  
  # Test finding best subset
  cat("\nBest subset (maximization):\n")
  best <- find_best_subset(subsets, dist_mat, maximize = TRUE)
  print(best)
  
  # Test minimization
  cat("\nRanking subsets (minimization):\n")
  ranked_min <- rank_subsets_lexico(subsets, dist_mat, maximize = FALSE)
  print(ranked_min)
  
  cat("\nBest subset (minimization):\n")
  best_min <- find_best_subset(subsets, dist_mat, maximize = FALSE)
  print(best_min)
}
