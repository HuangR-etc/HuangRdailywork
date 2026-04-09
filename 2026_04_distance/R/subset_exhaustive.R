# Exhaustive search for small trees
# This module implements exhaustive enumeration for validation on small trees

library(utils)  # For combn function

#' Generate all possible subsets of given size
#'
#' @param n_tips Total number of tips
#' @param subset_size Desired subset size
#' @return A matrix where each column is a subset (indices)
generate_all_subsets <- function(n_tips, subset_size) {
  if (subset_size > n_tips) {
    stop("Subset size cannot be larger than total number of tips")
  }
  
  # Generate all combinations
  all_combinations <- combn(n_tips, subset_size)
  
  return(all_combinations)
}

#' Find exact optimum using exhaustive search
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @param maximize If TRUE, find maximum; if FALSE, find minimum
#' @return A list containing the optimal subset and its metrics
find_exact_optimum <- function(dist_obj, subset_size, maximize = TRUE) {
  n_tips <- length(dist_obj$tip_labels)
  
  cat("Finding exact optimum for", n_tips, "tips, subset size", subset_size, "\n")
  cat("Total combinations: choose(", n_tips, ",", subset_size, ") =", 
      choose(n_tips, subset_size), "\n")
  
  # Generate all subsets
  all_subsets <- generate_all_subsets(n_tips, subset_size)
  n_combinations <- ncol(all_subsets)
  
  cat("Enumerating", n_combinations, "combinations...\n")
  
  # Initialize tracking
  best_subset <- NULL
  best_metrics <- NULL
  best_index <- 1
  
  # Progress tracking
  progress_interval <- max(1, floor(n_combinations / 10))
  
  # Evaluate all subsets
  for (i in 1:n_combinations) {
    # Progress reporting
    if (i %% progress_interval == 0) {
      cat("  Processed", i, "of", n_combinations, "combinations (", 
          round(i/n_combinations*100, 1), "%)\n")
    }
    
    # Get current subset
    current_subset <- all_subsets[, i]
    
    # Calculate metrics
    current_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
    
    # Check if this is the best so far
    if (is.null(best_subset)) {
      best_subset <- current_subset
      best_metrics <- current_metrics
      best_index <- i
    } else {
      if (maximize) {
        if (is_better_lexico_max(current_metrics, best_metrics)) {
          best_subset <- current_subset
          best_metrics <- current_metrics
          best_index <- i
        }
      } else {
        if (is_better_lexico_min(current_metrics, best_metrics)) {
          best_subset <- current_subset
          best_metrics <- current_metrics
          best_index <- i
        }
      }
    }
  }
  
  cat("Exact optimum found.\n")
  
  return(list(
    subset = best_subset,
    subset_names = dist_obj$tip_labels[best_subset],
    metrics = best_metrics,
    combination_index = best_index,
    total_combinations = n_combinations,
    algorithm = ifelse(maximize, "exact_max", "exact_min")
  ))
}

#' Compare heuristic result with exact optimum
#'
#' @param dist_obj Distance object
#' @param heuristic_result Result from run_complete_algorithm
#' @param subset_size Subset size
#' @param maximize If TRUE, compare for maximization (default: TRUE)
#' @return A comprehensive comparison
compare_heuristic_vs_exact <- function(dist_obj, heuristic_result, subset_size, maximize = TRUE) {
  cat("Comparing heuristic with exact optimum...\n")
  
  # Find exact optimum
  exact_result <- find_exact_optimum(dist_obj, subset_size, maximize)
  
  # Get heuristic metrics
  heuristic_metrics <- heuristic_result$final_metrics
  
  # Calculate gaps
  gaps <- list()
  gaps$MinPD <- ifelse(maximize, 
                       exact_result$metrics$MinPD - heuristic_metrics$MinPD,
                       heuristic_metrics$MinPD - exact_result$metrics$MinPD)
  gaps$MeanPD <- ifelse(maximize,
                        exact_result$metrics$MeanPD - heuristic_metrics$MeanPD,
                        heuristic_metrics$MeanPD - exact_result$metrics$MeanPD)
  gaps$MeanNND <- ifelse(maximize,
                         exact_result$metrics$MeanNND - heuristic_metrics$MeanNND,
                         heuristic_metrics$MeanNND - exact_result$metrics$MeanNND)
  
  # Calculate relative gaps (as percentage)
  rel_gaps <- list()
  rel_gaps$MinPD <- ifelse(exact_result$metrics$MinPD != 0,
                           gaps$MinPD / exact_result$metrics$MinPD * 100,
                           NA)
  rel_gaps$MeanPD <- ifelse(exact_result$metrics$MeanPD != 0,
                            gaps$MeanPD / exact_result$metrics$MeanPD * 100,
                            NA)
  rel_gaps$MeanNND <- ifelse(exact_result$metrics$MeanNND != 0,
                             gaps$MeanNND / exact_result$metrics$MeanNND * 100,
                             NA)
  
  # Check if heuristic found the exact optimum
  heuristic_subset_sorted <- sort(heuristic_result$final_subset)
  exact_subset_sorted <- sort(exact_result$subset)
  is_exact_match <- identical(heuristic_subset_sorted, exact_subset_sorted)
  
  # Create comparison summary
  comparison_summary <- data.frame(
    Metric = c("MinPD", "MeanPD", "MeanNND"),
    Heuristic = c(heuristic_metrics$MinPD, heuristic_metrics$MeanPD, heuristic_metrics$MeanNND),
    Exact = c(exact_result$metrics$MinPD, exact_result$metrics$MeanPD, exact_result$metrics$MeanNND),
    Gap = c(gaps$MinPD, gaps$MeanPD, gaps$MeanNND),
    Rel_Gap_Pct = c(rel_gaps$MinPD, rel_gaps$MeanPD, rel_gaps$MeanNND),
    stringsAsFactors = FALSE
  )
  
  return(list(
    heuristic_result = heuristic_result,
    exact_result = exact_result,
    comparison_summary = comparison_summary,
    gaps = gaps,
    rel_gaps = rel_gaps,
    is_exact_match = is_exact_match,
    heuristic_subset = heuristic_result$final_subset_names,
    exact_subset = exact_result$subset_names
  ))
}

#' Evaluate heuristic performance on small tree
#'
#' @param dist_obj Distance object (for small tree)
#' @param subset_size Subset size
#' @param maximize If TRUE, evaluate for maximization (default: TRUE)
#' @param n_greedy_starts Number of random starts for greedy phase
#' @return Comprehensive evaluation results
evaluate_heuristic_performance <- function(dist_obj, subset_size, maximize = TRUE, n_greedy_starts = 1) {
  cat("Evaluating heuristic performance on small tree...\n")
  
  # Run heuristic algorithm
  heuristic_result <- run_complete_algorithm(dist_obj, subset_size, maximize, 
                                            n_greedy_starts = n_greedy_starts)
  
  # Compare with exact optimum
  comparison <- compare_heuristic_vs_exact(dist_obj, heuristic_result, subset_size, maximize)
  
  # Calculate rank of heuristic solution among all possible subsets
  cat("Calculating rank of heuristic solution...\n")
  n_tips <- length(dist_obj$tip_labels)
  all_subsets <- generate_all_subsets(n_tips, subset_size)
  n_combinations <- ncol(all_subsets)
  
  # We'll sample if there are too many combinations
  max_evaluate <- 10000  # Limit for performance
  if (n_combinations > max_evaluate) {
    cat("  Too many combinations (", n_combinations, "), sampling", max_evaluate, "subsets\n")
    sampled_indices <- sample(1:n_combinations, max_evaluate)
    all_subsets <- all_subsets[, sampled_indices, drop = FALSE]
    n_combinations <- max_evaluate
  }
  
  # Count how many subsets are better than heuristic
  better_count <- 0
  equal_count <- 0
  
  for (i in 1:n_combinations) {
    current_subset <- all_subsets[, i]
    current_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
    
    if (maximize) {
      if (is_better_lexico_max(current_metrics, heuristic_result$final_metrics)) {
        better_count <- better_count + 1
      } else if (!is_better_lexico_max(heuristic_result$final_metrics, current_metrics) &&
                 !is_better_lexico_max(current_metrics, heuristic_result$final_metrics)) {
        equal_count <- equal_count + 1
      }
    } else {
      if (is_better_lexico_min(current_metrics, heuristic_result$final_metrics)) {
        better_count <- better_count + 1
      } else if (!is_better_lexico_min(heuristic_result$final_metrics, current_metrics) &&
                 !is_better_lexico_min(current_metrics, heuristic_result$final_metrics)) {
        equal_count <- equal_count + 1
      }
    }
  }
  
  # Calculate approximate rank and percentile
  rank <- better_count + 1
  percentile <- (n_combinations - better_count) / n_combinations * 100
  
  comparison$rank_analysis <- list(
    better_count = better_count,
    equal_count = equal_count,
    evaluated_count = n_combinations,
    rank = rank,
    percentile = percentile
  )
  
  return(comparison)
}

#' Test exhaustive search on a small tree
test_exhaustive <- function() {
  # Load required libraries and functions
  library(ape)
  source("distance_metrics.R")
  source("objective_compare.R")
  source("subset_greedy.R")
  source("subset_exchange.R")
  
  # Create a small test tree
  test_tree <- rtree(15)
  test_tree$tip.label <- paste0("sp", 1:15)
  
  # Create distance object
  dist_obj <- create_distance_object(test_tree)
  
  cat("Testing exhaustive search on small tree (15 tips):\n")
  
  # Test 1: Generate all subsets
  cat("\n1. Generating all subsets of size 5:\n")
  all_subsets <- generate_all_subsets(15, 5)
  cat("   Total combinations:", ncol(all_subsets), "\n")
  cat("   First subset:", all_subsets[, 1], "\n")
  cat("   First subset names:", dist_obj$tip_labels[all_subsets[, 1]], "\n")
  
  # Test 2: Find exact optimum (maximization)
  cat("\n2. Finding exact optimum (maximization):\n")
  exact_max <- find_exact_optimum(dist_obj, subset_size = 5, maximize = TRUE)
  cat("   Exact optimum subset:", exact_max$subset_names, "\n")
  cat("   Exact optimum metrics: MinPD =", exact_max$metrics$MinPD,
      "MeanPD =", exact_max$metrics$MeanPD,
      "MeanNND =", exact_max$metrics$MeanNND, "\n")
  
  # Test 3: Find exact optimum (minimization)
  cat("\n3. Finding exact optimum (minimization):\n")
  exact_min <- find_exact_optimum(dist_obj, subset_size = 5, maximize = FALSE)
  cat("   Exact optimum subset:", exact_min$subset_names, "\n")
  cat("   Exact optimum metrics: MinPD =", exact_min$metrics$MinPD,
      "MeanPD =", exact_min$metrics$MeanPD,
      "MeanNND =", exact_min$metrics$MeanNND, "\n")
  
  # Test 4: Run heuristic and compare
  cat("\n4. Running heuristic and comparing with exact optimum:\n")
  heuristic_result <- run_complete_algorithm(dist_obj, subset_size = 5, maximize = TRUE)
  comparison <- compare_heuristic_vs_exact(dist_obj, heuristic_result, subset_size = 5, maximize = TRUE)
  
  cat("   Heuristic subset:", comparison$heuristic_subset, "\n")
  cat("   Exact subset:", comparison$exact_subset, "\n")
  cat("   Exact match?", comparison$is_exact_match, "\n")
  cat("   Comparison summary:\n")
  print(comparison$comparison_summary)
  
  # Test 5: Full evaluation
  cat("\n5. Full heuristic performance evaluation:\n")
  evaluation <- evaluate_heuristic_performance(dist_obj, subset_size = 5, maximize = TRUE)
  cat("   Rank analysis:\n")
  cat("     Better subsets:", evaluation$rank_analysis$better_count, "\n")
  cat("     Equal subsets:", evaluation$rank_analysis$equal_count, "\n")
  cat("     Evaluated subsets:", evaluation$rank_analysis$evaluated_count, "\n")
  cat("     Rank:", evaluation$rank_analysis$rank, "\n")
  cat("     Percentile:", round(evaluation$rank_analysis$percentile, 2), "%\n")
  
  return(list(
    exact_max = exact_max,
    exact_min = exact_min,
    comparison = comparison,
    evaluation = evaluation
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_exhaustive()
}
