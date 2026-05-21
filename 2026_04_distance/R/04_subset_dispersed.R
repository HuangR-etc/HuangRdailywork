# 04_subset_dispersed.R
# Dispersed subset selection: fixed-start greedy + exchange refinement
# Migrated from old objective_compare.R, subset_greedy.R, subset_exchange.R

# ============================================================
# Lexicographic comparison (maximization only for dispersed)
# ============================================================

#' Compare two subsets using lexicographic ordering (maximization)
#'
#' @param metrics_a List with MinPD, MeanPD, MeanNND for subset A
#' @param metrics_b List with MinPD, MeanPD, MeanNND for subset B
#' @return TRUE if subset A is better than subset B, FALSE otherwise
is_better_lexico_max <- function(metrics_a, metrics_b) {
  if (metrics_a$MinPD > metrics_b$MinPD) return(TRUE)
  if (metrics_a$MinPD < metrics_b$MinPD) return(FALSE)
  if (metrics_a$MeanPD > metrics_b$MeanPD) return(TRUE)
  if (metrics_a$MeanPD < metrics_b$MeanPD) return(FALSE)
  if (metrics_a$MeanNND > metrics_b$MeanNND) return(TRUE)
  FALSE
}

# ============================================================
# Fixed-start dispersed greedy
# ============================================================

#' Find the species with the maximum mean distance to all others
#'
#' @param dist_obj Distance object
#' @return List with start_idx, start_name, mean_distance, all_mean_distances
find_max_mean_distance_species <- function(dist_obj) {
  d <- dist_obj$dist_mat
  diag(d) <- NA
  mean_dist <- rowMeans(d, na.rm = TRUE)
  start_idx <- which.max(mean_dist)
  list(
    start_idx = start_idx,
    start_name = dist_obj$tip_labels[start_idx],
    mean_distance = mean_dist[start_idx],
    all_mean_distances = mean_dist
  )
}

#' Build a dispersed subset using fixed-start greedy (maximize lexicographic)
#'
#' Starts from the species with highest mean distance to all others,
#' then greedily adds species that maximize MinPD > MeanPD > MeanNND.
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @return List with subset, subset_names, metrics, start_info, algorithm
build_dispersed_greedy_fixed_start <- function(dist_obj, subset_size) {
  all_tips <- seq_along(dist_obj$tip_labels)
  
  start_info <- find_max_mean_distance_species(dist_obj)
  current_subset <- start_info$start_idx
  available <- setdiff(all_tips, current_subset)
  
  while (length(current_subset) < subset_size) {
    best_candidate <- NULL
    best_metrics <- NULL
    
    for (candidate in available) {
      temp_subset <- c(current_subset, candidate)
      temp_metrics <- calc_subset_metrics(dist_obj$dist_mat, temp_subset)
      
      if (is.null(best_candidate)) {
        best_candidate <- candidate
        best_metrics <- temp_metrics
      } else if (is_better_lexico_max(temp_metrics, best_metrics)) {
        best_candidate <- candidate
        best_metrics <- temp_metrics
      }
    }
    
    current_subset <- c(current_subset, best_candidate)
    available <- setdiff(all_tips, current_subset)
  }
  
  final_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
  
  list(
    subset = current_subset,
    subset_names = dist_obj$tip_labels[current_subset],
    metrics = final_metrics,
    start_info = start_info,
    algorithm = "dispersed_fixed_max_mean_distance_greedy"
  )
}

# ============================================================
# Exchange refinement (maximization)
# ============================================================

#' Refine a subset using 1-for-1 exchange (maximization)
#'
#' @param dist_obj Distance object
#' @param current_subset Current subset (indices)
#' @param maximize If TRUE, maximize metrics; if FALSE, minimize metrics
#' @param single_objective If not NULL, optimize only this objective
#' @param max_iterations Maximum number of iterations
#' @return A list containing the refined subset and improvement information
refine_subset_exchange <- function(dist_obj, current_subset, maximize = TRUE,
                                   single_objective = NULL, max_iterations = 100) {
  all_tips <- seq_along(dist_obj$tip_labels)
  current_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
  
  improved <- TRUE
  iterations <- 0
  improvements <- list()
  
  while (improved && iterations < max_iterations) {
    improved <- FALSE
    selected <- current_subset
    available <- setdiff(all_tips, selected)
    
    for (out_species in selected) {
      for (in_species in available) {
        new_subset <- setdiff(selected, out_species)
        new_subset <- c(new_subset, in_species)
        
        if (is.null(single_objective)) {
          new_metrics <- calc_subset_metrics(dist_obj$dist_mat, new_subset)
          if (maximize) {
            if (is_better_lexico_max(new_metrics, current_metrics)) {
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              break
            }
          } else {
            # minimization not used for dispersed, but kept for completeness
            if (is_better_lexico_min(new_metrics, current_metrics)) {
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              break
            }
          }
        } else {
          new_metrics <- calc_subset_metrics(dist_obj$dist_mat, new_subset)
          new_value <- new_metrics[[single_objective]]
          current_value <- current_metrics[[single_objective]]
          if (maximize) {
            if (new_value > current_value) {
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              break
            }
          } else {
            if (new_value < current_value) {
              current_subset <- new_subset
              current_metrics <- new_metrics
              improved <- TRUE
              improvements[[length(improvements) + 1]] <- list(
                iteration = iterations + 1,
                out_species = out_species,
                in_species = in_species,
                metrics = new_metrics
              )
              break
            }
          }
        }
      }
      if (improved) break
    }
    iterations <- iterations + 1
  }
  
  list(
    subset = current_subset,
    subset_names = dist_obj$tip_labels[current_subset],
    metrics = current_metrics,
    improvements = improvements,
    iterations = iterations,
    converged = !improved,
    algorithm = ifelse(is.null(single_objective),
                      ifelse(maximize, "exchange_lexico_max", "exchange_lexico_min"),
                      paste0("exchange_", single_objective, ifelse(maximize, "_max", "_min")))
  )
}

# ============================================================
# Complete dispersed algorithm: fixed-start greedy + exchange
# ============================================================

#' Run the complete dispersed algorithm (fixed-start greedy + exchange)
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @return List with start_info, greedy_result, exchange_result, final_subset, etc.
run_dispersed_algorithm <- function(dist_obj, subset_size) {
  greedy_result <- build_dispersed_greedy_fixed_start(
    dist_obj = dist_obj,
    subset_size = subset_size
  )
  
  exchange_result <- refine_subset_exchange(
    dist_obj = dist_obj,
    current_subset = greedy_result$subset,
    maximize = TRUE,
    single_objective = NULL
  )
  
  list(
    start_info = greedy_result$start_info,
    greedy_result = greedy_result,
    exchange_result = exchange_result,
    final_subset = exchange_result$subset,
    final_subset_names = exchange_result$subset_names,
    final_metrics = calc_subset_metrics_extended(
      dist_obj$dist_mat,
      exchange_result$subset
    ),
    algorithm = "dispersed_fixed_start_greedy_exchange"
  )
}
