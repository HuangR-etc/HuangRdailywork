# 05_subset_clustered.R
# Clustered subset selection: multi-start greedy-plus-exchange method
# (primary method) and legacy seed-nearest-neighbor method (fast option)
#
# Primary method (new):
#   For each species as starting species:
#     1. Greedy forward selection minimizing MeanPD -> MeanNND -> MaxPD
#     2. One-for-one exchange refinement until no improvement
#   Then select best across all starting species by MeanPD -> MeanNND -> MaxPD
#
# Legacy method (old):
#   For each species as seed, take nearest (s - 1) neighbors.
#   Sort by MeanPD -> MeanNND -> MaxPD.

# ============================================================
# Lexicographic comparison for clustered minimization
# ============================================================

#' Compare two subsets under clustered lexicographic minimization
#'
#' Returns TRUE if metrics_a is better than metrics_b under:
#'   smaller MeanPD -> smaller MeanNND -> smaller MaxPD
#'
#' @param metrics_a List with MinPD, MeanPD, MeanNND, MaxPD
#' @param metrics_b List with MinPD, MeanPD, MeanNND, MaxPD (or NULL)
#' @param tol Tolerance for tie-breaking
#' @return TRUE if metrics_a is strictly better than metrics_b
is_better_lexico_min_clustered <- function(metrics_a, metrics_b, tol = 1e-10) {
  if (is.null(metrics_b)) return(TRUE)
  
  if (metrics_a$MeanPD < metrics_b$MeanPD - tol) return(TRUE)
  if (metrics_a$MeanPD > metrics_b$MeanPD + tol) return(FALSE)
  
  if (metrics_a$MeanNND < metrics_b$MeanNND - tol) return(TRUE)
  if (metrics_a$MeanNND > metrics_b$MeanNND + tol) return(FALSE)
  
  if (metrics_a$MaxPD < metrics_b$MaxPD - tol) return(TRUE)
  return(FALSE)
}

# ============================================================
# Greedy forward selection from a single starting species
# ============================================================

#' Build a clustered subset via greedy forward selection from one starting species
#'
#' Starting from start_idx, iteratively add the species that minimizes
#' MeanPD (then MeanNND, then MaxPD) until subset_size is reached.
#'
#' @param dist_obj Distance object (list with dist_mat, tip_labels)
#' @param subset_size Desired subset size s
#' @param start_idx Index of the starting species
#' @param tol Tolerance for tie-breaking
#' @return List with subset, subset_names, metrics, greedy_log, algorithm
build_clustered_greedy_from_start <- function(dist_obj,
                                              subset_size,
                                              start_idx,
                                              tol = 1e-10) {
  all_tips <- seq_along(dist_obj$tip_labels)
  selected <- start_idx
  available <- setdiff(all_tips, selected)
  
  greedy_log <- list()
  
  while (length(selected) < subset_size) {
    best_candidate <- NULL
    best_metrics <- NULL
    
    for (candidate in available) {
      candidate_subset <- c(selected, candidate)
      candidate_metrics <- calc_subset_metrics_extended(dist_obj$dist_mat, candidate_subset)
      
      if (is.null(best_candidate) ||
          is_better_lexico_min_clustered(candidate_metrics, best_metrics, tol)) {
        best_candidate <- candidate
        best_metrics <- candidate_metrics
      }
    }
    
    selected <- c(selected, best_candidate)
    available <- setdiff(available, best_candidate)
    
    greedy_log[[length(greedy_log) + 1]] <- list(
      step = length(selected),
      added_idx = best_candidate,
      added_name = dist_obj$tip_labels[best_candidate],
      metrics = best_metrics
    )
  }
  
  list(
    subset = selected,
    subset_names = dist_obj$tip_labels[selected],
    metrics = calc_subset_metrics_extended(dist_obj$dist_mat, selected),
    start_idx = start_idx,
    start_name = dist_obj$tip_labels[start_idx],
    greedy_log = greedy_log,
    algorithm = "clustered_greedy_from_start_meanpd_meannnd_maxpd"
  )
}

# ============================================================
# One-for-one exchange refinement
# ============================================================

#' Refine a clustered subset via one-for-one exchange
#'
#' Iteratively try all possible swaps of one selected species with one
#' unselected species. Accept the best-improving swap per iteration.
#' Stop when no improving swap can be found or max_iterations reached.
#'
#' @param dist_obj Distance object
#' @param current_subset Vector of indices of the current subset
#' @param max_iterations Maximum number of exchange iterations
#' @param tol Tolerance for tie-breaking
#' @return List with refined subset, metrics, exchange_log, iterations, converged
refine_clustered_exchange <- function(dist_obj,
                                      current_subset,
                                      max_iterations = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
                                      tol = CLUSTERED_EXCHANGE_TOL) {
  all_tips <- seq_along(dist_obj$tip_labels)
  current_subset <- unique(current_subset)
  current_metrics <- calc_subset_metrics_extended(dist_obj$dist_mat, current_subset)
  
  iteration <- 0
  exchange_log <- list()
  converged <- FALSE
  
  repeat {
    if (iteration >= max_iterations) break
    
    selected <- current_subset
    available <- setdiff(all_tips, selected)
    
    best_swap <- NULL
    best_subset <- NULL
    best_metrics <- NULL
    
    for (out_species in selected) {
      for (in_species in available) {
        candidate_subset <- c(setdiff(selected, out_species), in_species)
        candidate_metrics <- calc_subset_metrics_extended(dist_obj$dist_mat, candidate_subset)
        
        if (is_better_lexico_min_clustered(candidate_metrics, current_metrics, tol)) {
          if (is.null(best_metrics) ||
              is_better_lexico_min_clustered(candidate_metrics, best_metrics, tol)) {
            best_swap <- list(out_species = out_species, in_species = in_species)
            best_subset <- candidate_subset
            best_metrics <- candidate_metrics
          }
        }
      }
    }
    
    if (is.null(best_swap)) {
      converged <- TRUE
      break
    }
    
    iteration <- iteration + 1
    current_subset <- best_subset
    current_metrics <- best_metrics
    
    exchange_log[[iteration]] <- list(
      iteration = iteration,
      out_species = best_swap$out_species,
      in_species = best_swap$in_species,
      out_name = dist_obj$tip_labels[best_swap$out_species],
      in_name = dist_obj$tip_labels[best_swap$in_species],
      metrics = current_metrics
    )
  }
  
  list(
    subset = current_subset,
    subset_names = dist_obj$tip_labels[current_subset],
    metrics = current_metrics,
    exchange_log = exchange_log,
    iterations = iteration,
    converged = converged,
    algorithm = "clustered_exchange_lexico_min_meanpd_meannnd_maxpd"
  )
}

# ============================================================
# Multi-start greedy-plus-exchange main function
# ============================================================

#' Select clustered subset via multi-start greedy-plus-exchange
#'
#' For each species as starting species:
#'   1. Greedy forward selection minimizing MeanPD -> MeanNND -> MaxPD
#'   2. One-for-one exchange refinement until convergence
#' Then select the best final subset across all starts.
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size s
#' @param max_iterations Maximum exchange iterations per start
#' @param tol Tolerance for tie-breaking
#' @param verbose If TRUE, print progress messages
#' @return List with final_subset, final_metrics, candidate_metrics, etc.
select_clustered_greedy_exchange <- function(dist_obj,
                                             subset_size,
                                             max_iterations = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
                                             tol = CLUSTERED_EXCHANGE_TOL,
                                             verbose = FALSE) {
  all_tips <- seq_along(dist_obj$tip_labels)
  n_tips <- length(all_tips)
  
  if (subset_size >= n_tips) {
    stop("subset_size must be smaller than N.")
  }
  
  start_results <- vector("list", n_tips)
  
  for (start_idx in all_tips) {
    if (verbose) {
      cat("  Starting species", start_idx, "/", n_tips, ":",
          dist_obj$tip_labels[start_idx], "\n")
    }
    
    # Step 1: Greedy forward selection
    greedy_result <- build_clustered_greedy_from_start(
      dist_obj = dist_obj,
      subset_size = subset_size,
      start_idx = start_idx,
      tol = tol
    )
    
    # Step 2: Exchange refinement
    exchange_result <- refine_clustered_exchange(
      dist_obj = dist_obj,
      current_subset = greedy_result$subset,
      max_iterations = max_iterations,
      tol = tol
    )
    
    start_results[[start_idx]] <- list(
      start_idx = start_idx,
      start_name = dist_obj$tip_labels[start_idx],
      greedy_result = greedy_result,
      exchange_result = exchange_result,
      final_subset = exchange_result$subset,
      final_subset_names = exchange_result$subset_names,
      final_metrics = exchange_result$metrics
    )
  }
  
  # Build candidate metrics data frame
  metrics_df <- data.frame(
    Start_Idx = integer(n_tips),
    Start_Name = character(n_tips),
    Final_MinPD = numeric(n_tips),
    Final_MeanPD = numeric(n_tips),
    Final_MeanNND = numeric(n_tips),
    Final_MaxPD = numeric(n_tips),
    Greedy_MinPD = numeric(n_tips),
    Greedy_MeanPD = numeric(n_tips),
    Greedy_MeanNND = numeric(n_tips),
    Greedy_MaxPD = numeric(n_tips),
    Exchange_Iterations = integer(n_tips),
    Exchange_Converged = logical(n_tips),
    MinPD = numeric(n_tips),
    MeanPD = numeric(n_tips),
    MeanNND = numeric(n_tips),
    MaxPD = numeric(n_tips),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n_tips)) {
    sr <- start_results[[i]]
    metrics_df$Start_Idx[i] <- sr$start_idx
    metrics_df$Start_Name[i] <- sr$start_name
    metrics_df$Final_MinPD[i] <- sr$final_metrics$MinPD
    metrics_df$Final_MeanPD[i] <- sr$final_metrics$MeanPD
    metrics_df$Final_MeanNND[i] <- sr$final_metrics$MeanNND
    metrics_df$Final_MaxPD[i] <- sr$final_metrics$MaxPD
    metrics_df$Greedy_MinPD[i] <- sr$greedy_result$metrics$MinPD
    metrics_df$Greedy_MeanPD[i] <- sr$greedy_result$metrics$MeanPD
    metrics_df$Greedy_MeanNND[i] <- sr$greedy_result$metrics$MeanNND
    metrics_df$Greedy_MaxPD[i] <- sr$greedy_result$metrics$MaxPD
    metrics_df$Exchange_Iterations[i] <- sr$exchange_result$iterations
    metrics_df$Exchange_Converged[i] <- sr$exchange_result$converged
    # Compatibility columns for sorting
    metrics_df$MinPD[i] <- sr$final_metrics$MinPD
    metrics_df$MeanPD[i] <- sr$final_metrics$MeanPD
    metrics_df$MeanNND[i] <- sr$final_metrics$MeanNND
    metrics_df$MaxPD[i] <- sr$final_metrics$MaxPD
  }
  
  # Add Final_Subset_Key for deduplication
  metrics_df$Final_Subset_Key <- sapply(start_results, function(sr) {
    paste(sort(sr$final_subset), collapse = "|")
  })
  
  # Sort by MeanPD -> MeanNND -> MaxPD
  ord <- order(metrics_df$MeanPD, metrics_df$MeanNND, metrics_df$MaxPD)
  ordered_metrics <- metrics_df[ord, ]
  
  best_row <- metrics_df[ord[1], ]
  best_result <- start_results[[ord[1]]]
  
  # Count unique final subsets
  n_unique_final_subsets <- length(unique(metrics_df$Final_Subset_Key))
  
  list(
    final_subset = best_result$final_subset,
    final_subset_names = best_result$final_subset_names,
    final_metrics = best_result$final_metrics,
    
    candidate_metrics = metrics_df,
    ordered_candidate_metrics = ordered_metrics,
    
    start_results = start_results,
    n_raw_candidates = n_tips,
    n_unique_candidates = n_unique_final_subsets,
    
    best_seed_idx = best_row$Start_Idx,
    best_seed_name = best_row$Start_Name,
    best_start_idx = best_row$Start_Idx,
    best_start_name = best_row$Start_Name,
    
    greedy_result = best_result$greedy_result,
    exchange_result = best_result$exchange_result,
    
    algorithm = "clustered_multistart_greedy_exchange_meanpd_meannnd_maxpd"
  )
}

# ============================================================
# Legacy seed-nearest-neighbor method (fast option)
# ============================================================

#' Generate clustered candidate subsets via seed-neighborhood (legacy)
#'
#' For each species as seed, take the nearest (s - 1) neighbors.
#' Remove duplicate subsets (same set of species).
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @return List with subsets (unique), seed_idx, seed_names, n_raw, n_unique
generate_clustered_seed_neighborhoods <- function(dist_obj, subset_size) {
  d <- dist_obj$dist_mat
  all_tips <- seq_along(dist_obj$tip_labels)
  
  candidates <- vector("list", length(all_tips))
  
  for (seed in all_tips) {
    distances_from_seed <- d[seed, ]
    distances_from_seed[seed] <- Inf
    
    nearest_neighbors <- order(distances_from_seed)[1:(subset_size - 1)]
    subset_now <- c(seed, nearest_neighbors)
    
    candidates[[seed]] <- sort(unique(subset_now))
  }
  
  keys <- sapply(candidates, function(x) paste(sort(x), collapse = "|"))
  unique_keys <- unique(keys)
  unique_candidates <- candidates[match(unique_keys, keys)]
  seed_id_for_unique <- match(unique_keys, keys)
  
  list(
    subsets = unique_candidates,
    seed_idx = seed_id_for_unique,
    seed_names = dist_obj$tip_labels[seed_id_for_unique],
    n_raw = length(candidates),
    n_unique = length(unique_candidates)
  )
}

#' Select the best clustered neighborhood (legacy nearest-neighbor method)
#'
#' Sort candidates by MeanPD (asc), MeanNND (asc), MaxPD (asc).
#' Return the best (first) candidate.
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @return List with final_subset, final_subset_names, final_metrics,
#'         candidate_metrics, ordered_candidate_metrics, best_seed info
select_clustered_nearest_neighbor <- function(dist_obj, subset_size) {
  cand <- generate_clustered_seed_neighborhoods(
    dist_obj = dist_obj,
    subset_size = subset_size
  )
  
  metrics_df <- calc_multiple_subsets_metrics_extended(
    dist_mat = dist_obj$dist_mat,
    subsets = cand$subsets
  )
  
  metrics_df$Seed_Idx <- cand$seed_idx
  metrics_df$Seed_Name <- cand$seed_names
  
  # Sort by MeanPD (asc), MeanNND (asc), MaxPD (asc)
  ord <- order(
    metrics_df$MeanPD,
    metrics_df$MeanNND,
    metrics_df$MaxPD
  )
  
  best_row <- metrics_df[ord[1], ]
  best_subset <- cand$subsets[[ord[1]]]
  
  list(
    final_subset = best_subset,
    final_subset_names = dist_obj$tip_labels[best_subset],
    final_metrics = calc_subset_metrics_extended(dist_obj$dist_mat, best_subset),
    candidate_metrics = metrics_df,
    ordered_candidate_metrics = metrics_df[ord, ],
    n_raw_candidates = cand$n_raw,
    n_unique_candidates = cand$n_unique,
    best_seed_idx = best_row$Seed_Idx,
    best_seed_name = best_row$Seed_Name,
    algorithm = "clustered_seed_nearest_neighbors_meanpd_meannnd_maxpd"
  )
}

# ============================================================
# Deprecated wrapper (backward compatibility)
# ============================================================

#' Deprecated: use select_clustered_greedy_exchange() for main analyses
#'
#' This wrapper calls the legacy nearest-neighbor method.
#' It is kept for backward compatibility but should not be used
#' for main analyses.
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @param ... Additional arguments (ignored)
#' @return Result from select_clustered_nearest_neighbor()
select_best_clustered_neighborhood <- function(dist_obj, subset_size, ...) {
  warning("select_best_clustered_neighborhood() is deprecated; ",
          "use select_clustered_greedy_exchange() for main analyses.")
  select_clustered_nearest_neighbor(dist_obj, subset_size)
}
