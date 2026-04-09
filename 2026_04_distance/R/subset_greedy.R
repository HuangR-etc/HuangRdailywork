# Greedy construction algorithm for phylogenetic dispersed subset analysis
# This module implements the greedy construction phase of the main algorithm

#' Build a subset using greedy construction
#'
#' @param dist_obj Distance object (from create_distance_object)
#' @param subset_size Desired subset size
#' @param maximize If TRUE, maximize metrics; if FALSE, minimize metrics
#' @param single_objective If not NULL, optimize only this objective ("MinPD", "MeanPD", or "MeanNND")
#' @param start_species Optional starting species (index or name). If NULL, random start.
#' @return A list containing the subset and metrics
build_subset_greedy <- function(dist_obj, subset_size, maximize = TRUE, 
                                single_objective = NULL, start_species = NULL) {
  
  # Get all tip indices
  all_tips <- seq_along(dist_obj$tip_labels)
  
  # Initialize subset
  if (is.null(start_species)) {
    # Random start
    current_subset <- sample(all_tips, 1)
  } else {
    # Use provided start species
    if (is.character(start_species)) {
      current_subset <- which(dist_obj$tip_labels %in% start_species)
    } else {
      current_subset <- start_species
    }
    if (length(current_subset) == 0) {
      stop("Start species not found in tree")
    }
  }
  
  # Get available species (not in current subset)
  available <- setdiff(all_tips, current_subset)
  
  # Greedy construction loop
  while (length(current_subset) < subset_size && length(available) > 0) {
    best_candidate <- NULL
    best_metrics <- NULL
    
    # Try each available candidate
    for (candidate in available) {
      # Create temporary subset with candidate added
      temp_subset <- c(current_subset, candidate)
      
      # Calculate metrics for temporary subset
      if (is.null(single_objective)) {
        # Use lexicographic ordering
        temp_metrics <- calc_subset_metrics(dist_obj$dist_mat, temp_subset)
        
        # Compare with current best
        if (is.null(best_candidate)) {
          best_candidate <- candidate
          best_metrics <- temp_metrics
        } else {
          if (maximize) {
            if (is_better_lexico_max(temp_metrics, best_metrics)) {
              best_candidate <- candidate
              best_metrics <- temp_metrics
            }
          } else {
            if (is_better_lexico_min(temp_metrics, best_metrics)) {
              best_candidate <- candidate
              best_metrics <- temp_metrics
            }
          }
        }
      } else {
        # Single objective optimization
        temp_metrics <- calc_subset_metrics(dist_obj$dist_mat, temp_subset)
        temp_value <- temp_metrics[[single_objective]]
        
        if (is.null(best_candidate)) {
          best_candidate <- candidate
          best_value <- temp_value
        } else {
          if (maximize) {
            if (temp_value > best_value) {
              best_candidate <- candidate
              best_value <- temp_value
            }
          } else {
            if (temp_value < best_value) {
              best_candidate <- candidate
              best_value <- temp_value
            }
          }
        }
      }
    }
    
    # Add best candidate to subset
    current_subset <- c(current_subset, best_candidate)
    available <- setdiff(available, best_candidate)
  }
  
  # Calculate final metrics
  final_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
  
  # Return result
  return(list(
    subset = current_subset,
    subset_names = dist_obj$tip_labels[current_subset],
    metrics = final_metrics,
    algorithm = ifelse(is.null(single_objective), 
                      ifelse(maximize, "greedy_lexico_max", "greedy_lexico_min"),
                      paste0("greedy_", single_objective, ifelse(maximize, "_max", "_min")))
  ))
}

#' Run greedy construction with multiple random starts
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @param n_starts Number of random starts
#' @param maximize If TRUE, maximize metrics; if FALSE, minimize metrics
#' @param single_objective If not NULL, optimize only this objective
#' @return The best subset across all random starts
build_subset_greedy_multistart <- function(dist_obj, subset_size, n_starts = 10,
                                          maximize = TRUE, single_objective = NULL) {
  
  best_result <- NULL
  best_metrics <- NULL
  
  for (i in 1:n_starts) {
    # Run greedy construction with random start
    result <- build_subset_greedy(dist_obj, subset_size, maximize, 
                                  single_objective, start_species = NULL)
    
    # Compare with current best
    if (is.null(best_result)) {
      best_result <- result
      best_metrics <- result$metrics
    } else {
      if (maximize) {
        if (is_better_lexico_max(result$metrics, best_metrics)) {
          best_result <- result
          best_metrics <- result$metrics
        }
      } else {
        if (is_better_lexico_min(result$metrics, best_metrics)) {
          best_result <- result
          best_metrics <- result$metrics
        }
      }
    }
  }
  
  return(best_result)
}

#' Test greedy construction on a simple example
test_greedy <- function() {
  # Load required libraries and functions
  library(ape)
  source("distance_metrics.R")
  source("objective_compare.R")
  
  # Create a test tree
  test_tree <- rtree(20)
  test_tree$tip.label <- paste0("sp", 1:20)
  
  # Create distance object
  dist_obj <- create_distance_object(test_tree)
  
  cat("Testing greedy construction (maximization, lexicographic):\n")
  result_max <- build_subset_greedy(dist_obj, subset_size = 5, maximize = TRUE)
  cat("  Subset size:", length(result_max$subset), "\n")
  cat("  Subset:", result_max$subset_names, "\n")
  cat("  Metrics: MinPD =", result_max$metrics$MinPD, 
      "MeanPD =", result_max$metrics$MeanPD, 
      "MeanNND =", result_max$metrics$MeanNND, "\n")
  
  cat("\nTesting greedy construction (minimization, lexicographic):\n")
  result_min <- build_subset_greedy(dist_obj, subset_size = 5, maximize = FALSE)
  cat("  Subset size:", length(result_min$subset), "\n")
  cat("  Subset:", result_min$subset_names, "\n")
  cat("  Metrics: MinPD =", result_min$metrics$MinPD, 
      "MeanPD =", result_min$metrics$MeanPD, 
      "MeanNND =", result_min$metrics$MeanNND, "\n")
  
  cat("\nTesting greedy construction (single objective, maximize MeanPD):\n")
  result_single <- build_subset_greedy(dist_obj, subset_size = 5, maximize = TRUE,
                                       single_objective = "MeanPD")
  cat("  Subset size:", length(result_single$subset), "\n")
  cat("  Subset:", result_single$subset_names, "\n")
  cat("  Metrics: MinPD =", result_single$metrics$MinPD, 
      "MeanPD =", result_single$metrics$MeanPD, 
      "MeanNND =", result_single$metrics$MeanNND, "\n")
  
  cat("\nTesting multi-start greedy (maximization, 5 starts):\n")
  result_multi <- build_subset_greedy_multistart(dist_obj, subset_size = 5, 
                                                 n_starts = 5, maximize = TRUE)
  cat("  Subset size:", length(result_multi$subset), "\n")
  cat("  Subset:", result_multi$subset_names, "\n")
  cat("  Metrics: MinPD =", result_multi$metrics$MinPD, 
      "MeanPD =", result_multi$metrics$MeanPD, 
      "MeanNND =", result_multi$metrics$MeanNND, "\n")
  
  return(list(
    max_result = result_max,
    min_result = result_min,
    single_result = result_single,
    multi_result = result_multi
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_greedy()
}
