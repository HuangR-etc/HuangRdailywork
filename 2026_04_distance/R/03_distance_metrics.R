# 03_distance_metrics.R
# Distance metric calculations: MinPD, MeanPD, MeanNND, MaxPD
# Migrated from old distance_metrics.R, with MaxPD added.

#' Calculate patristic distance matrix for a tree
#'
#' @param tree A phylo object
#' @return A distance matrix
calc_distance_matrix <- function(tree) {
  ape::cophenetic.phylo(tree)
}

#' Calculate MinPD (minimum pairwise distance) for a subset
#'
#' @param dist_mat Distance matrix
#' @param subset Vector of tip indices or names
#' @return MinPD value
calc_minpd <- function(dist_mat, subset) {
  if (length(subset) < 2) {
    return(0)
  }
  sub_dist <- dist_mat[subset, subset]
  diag(sub_dist) <- NA
  min(sub_dist, na.rm = TRUE)
}

#' Calculate MeanPD (mean pairwise distance) for a subset
#'
#' @param dist_mat Distance matrix
#' @param subset Vector of tip indices or names
#' @return MeanPD value
calc_meanpd <- function(dist_mat, subset) {
  n <- length(subset)
  if (n < 2) {
    return(0)
  }
  sub_dist <- dist_mat[subset, subset]
  upper_tri <- sub_dist[upper.tri(sub_dist)]
  mean(upper_tri)
}

#' Calculate MeanNND (mean nearest neighbor distance) for a subset
#'
#' @param dist_mat Distance matrix
#' @param subset Vector of tip indices or names
#' @return MeanNND value
calc_meannnd <- function(dist_mat, subset) {
  n <- length(subset)
  if (n < 2) {
    return(0)
  }
  sub_dist <- dist_mat[subset, subset]
  diag(sub_dist) <- Inf
  nnd_values <- apply(sub_dist, 1, min)
  mean(nnd_values)
}

#' Calculate MaxPD (maximum pairwise distance) for a subset
#'
#' @param dist_mat Distance matrix
#' @param subset Vector of tip indices or names
#' @return MaxPD value
calc_maxpd <- function(dist_mat, subset) {
  if (length(subset) < 2) return(0)
  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  upper_values <- sub_dist[upper.tri(sub_dist)]
  max(upper_values, na.rm = TRUE)
}

#' Calculate all three base metrics for a subset (MinPD, MeanPD, MeanNND)
#'
#' @param dist_mat Distance matrix
#' @param subset Vector of tip indices or names
#' @return A list with MinPD, MeanPD, and MeanNND
calc_subset_metrics <- function(dist_mat, subset) {
  list(
    MinPD = calc_minpd(dist_mat, subset),
    MeanPD = calc_meanpd(dist_mat, subset),
    MeanNND = calc_meannnd(dist_mat, subset)
  )
}

#' Calculate extended metrics for a subset (MinPD, MeanPD, MeanNND, MaxPD)
#'
#' @param dist_mat Distance matrix
#' @param subset Vector of tip indices or names
#' @return A list with MinPD, MeanPD, MeanNND, MaxPD
calc_subset_metrics_extended <- function(dist_mat, subset) {
  base <- calc_subset_metrics(dist_mat, subset)
  list(
    MinPD = base$MinPD,
    MeanPD = base$MeanPD,
    MeanNND = base$MeanNND,
    MaxPD = calc_maxpd(dist_mat, subset)
  )
}

#' Calculate metrics for multiple subsets (base: MinPD, MeanPD, MeanNND)
#'
#' @param dist_mat Distance matrix
#' @param subsets List of subsets (each is a vector of tip indices/names)
#' @return A data frame with metrics for each subset
calc_multiple_subsets_metrics <- function(dist_mat, subsets) {
  results <- data.frame(
    SubsetID = integer(),
    MinPD = numeric(),
    MeanPD = numeric(),
    MeanNND = numeric(),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(subsets)) {
    metrics <- calc_subset_metrics(dist_mat, subsets[[i]])
    results <- rbind(results, data.frame(
      SubsetID = i,
      MinPD = metrics$MinPD,
      MeanPD = metrics$MeanPD,
      MeanNND = metrics$MeanNND,
      stringsAsFactors = FALSE
    ))
  }
  results
}

#' Calculate extended metrics for multiple subsets (MinPD, MeanPD, MeanNND, MaxPD)
#'
#' @param dist_mat Distance matrix
#' @param subsets List of subsets (each is a vector of tip indices/names)
#' @return A data frame with metrics for each subset
calc_multiple_subsets_metrics_extended <- function(dist_mat, subsets) {
  out <- data.frame(
    SubsetID = seq_along(subsets),
    MinPD = NA_real_,
    MeanPD = NA_real_,
    MeanNND = NA_real_,
    MaxPD = NA_real_
  )
  for (i in seq_along(subsets)) {
    m <- calc_subset_metrics_extended(dist_mat, subsets[[i]])
    out$MinPD[i] <- m$MinPD
    out$MeanPD[i] <- m$MeanPD
    out$MeanNND[i] <- m$MeanNND
    out$MaxPD[i] <- m$MaxPD
  }
  out
}

#' Create a distance object for efficient calculations
#'
#' @param tree A phylo object
#' @return A list containing tree and distance matrix
create_distance_object <- function(tree) {
  list(
    tree = tree,
    dist_mat = calc_distance_matrix(tree),
    tip_labels = tree$tip.label
  )
}

#' Get subset indices from tip names
#'
#' @param dist_obj Distance object
#' @param tip_names Vector of tip names
#' @return Vector of indices
get_subset_indices <- function(dist_obj, tip_names) {
  which(dist_obj$tip_labels %in% tip_names)
}

#' Get subset names from indices
#'
#' @param dist_obj Distance object
#' @param indices Vector of indices
#' @return Vector of tip names
get_subset_names <- function(dist_obj, indices) {
  dist_obj$tip_labels[indices]
}
