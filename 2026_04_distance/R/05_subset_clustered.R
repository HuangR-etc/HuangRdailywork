# 05_subset_clustered.R
# Clustered subset selection: seed-nearest-neighbor method
# NOT the old minimization greedy + exchange.

#' Generate clustered candidate subsets via seed-neighborhood
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

#' Select the best clustered neighborhood
#'
#' Sort candidates by MeanPD (asc), MeanNND (asc), MaxPD (asc).
#' Return the best (first) candidate.
#'
#' @param dist_obj Distance object
#' @param subset_size Desired subset size
#' @return List with final_subset, final_subset_names, final_metrics,
#'         candidate_metrics, ordered_candidate_metrics, best_seed info
select_best_clustered_neighborhood <- function(dist_obj, subset_size) {
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
