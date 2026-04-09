# Exchange refinement algorithm for phylogenetic dispersed subset analysis
# This module implements the exchange refinement phase of the main algorithm

#' Refine a subset using 1-for-1 exchange
#'
#' @param dist_obj Distance object
#' @param current_subset Current subset (indices)
#' @param maximize If TRUE, maximize metrics; if FALSE, minimize metrics
#' @param single_objective If not NULL, optimize only this objective
#' @param max_iterations Maximum number of iterations (default: 100)
#' @return A list containing the refined subset and improvement information
refine_subset_exchange <- function(dist_obj, current_subset, maximize = TRUE,
                                   single_objective = NULL, max_iterations = 100) {
  
  # Get all tip indices
  all_tips <- seq_along(dist_obj$tip_labels)
  
  # Get current metrics
  current_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
  
  # Initialize tracking
  improved <- TRUE
  iterations <- 0
  improvements <- list()
  
  while (improved && iterations < max_iterations) {
    improved <- FALSE
    
    # Get current selected and available species
    selected <- current_subset
    available <- setdiff(all_tips, selected)
    
    # Try all possible 1-for-1 exchanges
    for (out_species in selected) {
      for (in_species in available) {
        # Create new subset by swapping
        new_subset <- setdiff(selected, out_species)
        new_subset <- c(new_subset, in_species)
        
        # Calculate metrics for new subset
        if (is.null(single_objective)) {
          # Use lexicographic ordering
          new_metrics <- calc_subset_metrics(dist_obj$dist_mat, new_subset)
          
          # Check if new subset is better
          if (maximize) {
            if (is_better_lexico_max(new_metrics, current_metrics)) {
              # Accept the exchange
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              
              # Record improvement
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              
              # Break out of inner loops to restart scanning
              break
            }
          } else {
            if (is_better_lexico_min(new_metrics, current_metrics)) {
              # Accept the exchange
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              
              # Record improvement
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              
              # Break out of inner loops to restart scanning
              break
            }
          }
        } else {
          # Single objective optimization
          new_metrics <- calc_subset_metrics(dist_obj$dist_mat, new_subset)
          new_value <- new_metrics[[single_objective]]
          current_value <- current_metrics[[single_objective]]
          
          if (maximize) {
            if (new_value > current_value) {
              # Accept the exchange
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              
              # Record improvement
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              
              # Break out of inner loops to restart scanning
              break
            }
          } else {
            if (new_value < current_value) {
              # Accept the exchange
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              
              # Record improvement
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              
              # Break out of inner loops to restart scanning
              break
            }
          }
        }
      }
      
      if (improved) {
        break  # Break out of outer loop to restart scanning
      }
    }
    
    iterations <- iterations + 1
  }
  
  # Return result
  return(list(
    subset = current_subset,
    subset_names = dist_obj$tip_labels[current_subset],
    metrics = current_metrics,
    improvements = improvements,
    iterations = iterations,
    converged = !improved,
    algorithm = ifelse(is.null(single_objective), 
                      ifelse(maximize, "exchange_lexico_max", "exchange_lexico_min"),
                      paste0("exchange_", single_objective, ifelse(maximize, "_max", "_min")))
  ))
}

#' Run the complete algorithm (greedy + exchange)
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @param maximize If TRUE, maximize metrics; if FALSE, minimize metrics
#' @param single_objective If not NULL, optimize only this objective
#' @param n_greedy_starts Number of random starts for greedy phase
#' @return A list containing the final subset and algorithm details
run_complete_algorithm <- function(dist_obj, subset_size, maximize = TRUE,
                                   single_objective = NULL, n_greedy_starts = 1) {
  
  cat("Running complete algorithm...\n")
  
  # Phase 1: Greedy construction
  cat("  Phase 1: Greedy construction\n")
  if (n_greedy_starts > 1) {
    greedy_result <- build_subset_greedy_multistart(dist_obj, subset_size, 
                                                   n_starts = n_greedy_starts,
                                                   maximize = maximize,
                                                   single_objective = single_objective)
  } else {
    greedy_result <- build_subset_greedy(dist_obj, subset_size, 
                                        maximize = maximize,
                                        single_objective = single_objective)
  }
  
  cat("    Greedy result: MinPD =", greedy_result$metrics$MinPD,
      "MeanPD =", greedy_result$metrics$MeanPD,
      "MeanNND =", greedy_result$metrics$MeanNND, "\n")
  
  # Phase 2: Exchange refinement
  cat("  Phase 2: Exchange refinement\n")
  exchange_result <- refine_subset_exchange(dist_obj, greedy_result$subset,
                                           maximize = maximize,
                                           single_objective = single_objective)
  
  cat("    Exchange result: MinPD =", exchange_result$metrics$MinPD,
      "MeanPD =", exchange_result$metrics$MeanPD,
      "MeanNND =", exchange_result$metrics$MeanNND, "\n")
  cat("    Iterations:", exchange_result$iterations, "\n")
  cat("    Improvements:", length(exchange_result$improvements), "\n")
  
  # Calculate improvement from greedy to exchange
  improvement <- list()
  if (maximize) {
    improvement$MinPD <- exchange_result$metrics$MinPD - greedy_result$metrics$MinPD
    improvement$MeanPD <- exchange_result$metrics$MeanPD - greedy_result$metrics$MeanPD
    improvement$MeanNND <- exchange_result$metrics$MeanNND - greedy_result$metrics$MeanNND
  } else {
    improvement$MinPD <- greedy_result$metrics$MinPD - exchange_result$metrics$MinPD
    improvement$MeanPD <- greedy_result$metrics$MeanPD - exchange_result$metrics$MeanPD
    improvement$MeanNND <- greedy_result$metrics$MeanNND - exchange_result$metrics$MeanNND
  }
  
  # Return combined result
  return(list(
    greedy_result = greedy_result,
    exchange_result = exchange_result,
    final_subset = exchange_result$subset,
    final_subset_names = exchange_result$subset_names,
    final_metrics = exchange_result$metrics,
    improvement = improvement,
    algorithm = exchange_result$algorithm
  ))
}

#' Test exchange refinement on a simple example
test_exchange <- function() {
  # Load required libraries and functions
  library(ape)
  source("distance_metrics.R")
  source("objective_compare.R")
  source("subset_greedy.R")
  
  # Create a test tree
  test_tree <- rtree(20)
  test_tree$tip.label <- paste0("sp", 1:20)
  
  # Create distance object
  dist_obj <- create_distance_object(test_tree)
  
  # Create a starting subset (using greedy)
  cat("Testing exchange refinement:\n")
  
  # Test 1: Maximization with lexicographic ordering
  cat("\n1. Maximization with lexicographic ordering:\n")
  greedy_result <- build_subset_greedy(dist_obj, subset_size = 5, maximize = TRUE)
  cat("   Greedy subset:", greedy_result$subset_names, "\n")
  cat("   Greedy metrics: MinPD =", greedy_result$metrics$MinPD,
      "MeanPD =", greedy_result$metrics$MeanPD,
      "MeanNND =", greedy_result$metrics$MeanNND, "\n")
  
  exchange_result <- refine_subset_exchange(dist_obj, greedy_result$subset, maximize = TRUE)
  cat("   Exchange subset:", exchange_result$subset_names, "\n")
  cat("   Exchange metrics: MinPD =", exchange_result$metrics$MinPD,
      "MeanPD =", exchange_result$metrics$MeanPD,
      "MeanNND =", exchange_result$metrics$MeanNND, "\n")
  cat("   Iterations:", exchange_result$iterations, "\n")
  cat("   Improvements:", length(exchange_result$improvements), "\n")
  
  # Test 2: Complete algorithm
  cat("\n2. Complete algorithm (greedy + exchange):\n")
  complete_result <- run_complete_algorithm(dist_obj, subset_size = 5, maximize = TRUE)
  cat("   Final subset:", complete_result$final_subset_names, "\n")
  cat("   Final metrics: MinPD =", complete_result$final_metrics$MinPD,
      "MeanPD =", complete_result$final_metrics$MeanPD,
      "MeanNND =", complete_result$final_metrics$MeanNND, "\n")
  
  # Test 3: Minimization (for clustered subsets)
  cat("\n3. Minimization for clustered subsets:\n")
  complete_result_min <- run_complete_algorithm(dist_obj, subset_size = 5, maximize = FALSE)
  cat("   Final subset:", complete_result_min$final_subset_names, "\n")
  cat("   Final metrics: MinPD =", complete_result_min$final_metrics$MinPD,
      "MeanPD =", complete_result_min$final_metrics$MeanPD,
      "MeanNND =", complete_result_min$final_metrics$MeanNND, "\n")
  
  return(list(
    exchange_test = exchange_result,
    complete_test = complete_result,
    complete_test_min = complete_result_min
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_exchange()
}
