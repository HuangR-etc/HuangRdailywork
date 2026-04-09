# Distance metric calculations for phylogenetic dispersed subset analysis
# This module calculates MinPD, MeanPD, and MeanNND metrics

#' Calculate patristic distance matrix for a tree
#'
#' @param tree A phylo object
#' @return A distance matrix
calc_distance_matrix <- function(tree) {
  # Use cophenetic.phylo from ape package
  return(cophenetic.phylo(tree))
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
  
  # Get submatrix for the subset
  sub_dist <- dist_mat[subset, subset]
  
  # Set diagonal to NA to exclude self-distances
  diag(sub_dist) <- NA
  
  # Return minimum distance (excluding NA)
  return(min(sub_dist, na.rm = TRUE))
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
  
  # Get submatrix for the subset
  sub_dist <- dist_mat[subset, subset]
  
  # Get upper triangle (excluding diagonal)
  upper_tri <- sub_dist[upper.tri(sub_dist)]
  
  # Return mean
  return(mean(upper_tri))
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
  
  # Get submatrix for the subset
  sub_dist <- dist_mat[subset, subset]
  
  # Set diagonal to Inf to exclude self-distances
  diag(sub_dist) <- Inf
  
  # Calculate nearest neighbor distance for each species
  nnd_values <- apply(sub_dist, 1, min)
  
  # Return mean
  return(mean(nnd_values))
}

#' Calculate all three metrics for a subset
#'
#' @param dist_mat Distance matrix
#' @param subset Vector of tip indices or names
#' @return A list with MinPD, MeanPD, and MeanNND
calc_subset_metrics <- function(dist_mat, subset) {
  return(list(
    MinPD = calc_minpd(dist_mat, subset),
    MeanPD = calc_meanpd(dist_mat, subset),
    MeanNND = calc_meannnd(dist_mat, subset)
  ))
}

#' Calculate metrics for multiple subsets
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
  
  return(results)
}

#' Create a distance object for efficient calculations
#'
#' @param tree A phylo object
#' @return A list containing tree and distance matrix
create_distance_object <- function(tree) {
  return(list(
    tree = tree,
    dist_mat = calc_distance_matrix(tree),
    tip_labels = tree$tip.label
  ))
}

#' Get subset indices from tip names
#'
#' @param dist_obj Distance object
#' @param tip_names Vector of tip names
#' @return Vector of indices
get_subset_indices <- function(dist_obj, tip_names) {
  return(which(dist_obj$tip_labels %in% tip_names))
}

#' Get subset names from indices
#'
#' @param dist_obj Distance object
#' @param indices Vector of indices
#' @return Vector of tip names
get_subset_names <- function(dist_obj, indices) {
  return(dist_obj$tip_labels[indices])
}

# Test function
if (sys.nframe() == 0) {
  # Load required libraries
  library(ape)
  
  # Create a test tree
  test_tree <- rtree(10)
  test_tree$tip.label <- paste0("sp", 1:10)
  
  # Create distance object
  dist_obj <- create_distance_object(test_tree)
  
  # Test subset
  test_subset <- c("sp1", "sp3", "sp7", "sp9")
  test_indices <- get_subset_indices(dist_obj, test_subset)
  
  # Calculate metrics
  metrics <- calc_subset_metrics(dist_obj$dist_mat, test_indices)
  
  cat("Test tree with", length(test_tree$tip.label), "tips\n")
  cat("Test subset:", test_subset, "\n")
  cat("Metrics:\n")
  cat("  MinPD:", metrics$MinPD, "\n")
  cat("  MeanPD:", metrics$MeanPD, "\n")
  cat("  MeanNND:", metrics$MeanNND, "\n")
  
  # Test multiple subsets
  subsets <- list(
    c("sp1", "sp3", "sp5"),
    c("sp2", "sp4", "sp6"),
    c("sp1", "sp5", "sp9")
  )
  
  multi_metrics <- calc_multiple_subsets_metrics(dist_obj$dist_mat, 
    lapply(subsets, function(x) get_subset_indices(dist_obj, x)))
  
  cat("\nMultiple subsets metrics:\n")
  print(multi_metrics)
}
