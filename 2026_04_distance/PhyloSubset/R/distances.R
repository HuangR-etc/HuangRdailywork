#' Phylogenetic distance metrics for species subsets
#'
#' Functions to compute MinPD, MaxPD, MeanPD, MeanNND, and nearest neighbor
#' distances for a subset of species given a patristic distance matrix.
#'
#' @name distances
NULL

#' @rdname distances
#' @export
#' @title Compute all distance metrics for a subset
#' @description Calculates MinPD, MaxPD, MeanPD, and MeanNND for a given
#'   subset of species from a patristic distance matrix.
#' @param dist_mat A symmetric patristic distance matrix with species names
#'   as row and column names.
#' @param subset Character vector of species names in the subset.
#' @return A list with components:
#'   \item{MinPD}{Minimum pairwise distance.}
#'   \item{MaxPD}{Maximum pairwise distance.}
#'   \item{MeanPD}{Mean pairwise distance.}
#'   \item{MeanNND}{Mean nearest-neighbor distance.}
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(10)
#' D <- ape::cophenetic.phylo(tree)
#' subset <- sample(tree$tip.label, 4)
#' distance_metrics(D, subset)
#' }
distance_metrics <- function(dist_mat, subset) {
  if (length(subset) < 2) {
    return(list(MinPD = 0, MaxPD = 0, MeanPD = 0, MeanNND = 0))
  }

  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  diag(sub_dist) <- NA
  upper_vals <- sub_dist[upper.tri(sub_dist)]

  nnd <- apply(sub_dist, 1, min, na.rm = TRUE)

  list(
    MinPD = min(upper_vals, na.rm = TRUE),
    MaxPD = max(upper_vals, na.rm = TRUE),
    MeanPD = mean(upper_vals, na.rm = TRUE),
    MeanNND = mean(nnd, na.rm = TRUE)
  )
}

#' @rdname distances
#' @export
#' @title Minimum pairwise distance
#' @param dist_mat A symmetric patristic distance matrix.
#' @param subset Character vector of species names in the subset.
#' @return The minimum pairwise patristic distance.
min_pd <- function(dist_mat, subset) {
  if (length(subset) < 2) return(0)
  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  diag(sub_dist) <- NA
  min(sub_dist, na.rm = TRUE)
}

#' @rdname distances
#' @export
#' @title Maximum pairwise distance
#' @param dist_mat A symmetric patristic distance matrix.
#' @param subset Character vector of species names in the subset.
#' @return The maximum pairwise patristic distance.
max_pd <- function(dist_mat, subset) {
  if (length(subset) < 2) return(0)
  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  upper_vals <- sub_dist[upper.tri(sub_dist)]
  max(upper_vals, na.rm = TRUE)
}

#' @rdname distances
#' @export
#' @title Mean pairwise distance
#' @param dist_mat A symmetric patristic distance matrix.
#' @param subset Character vector of species names in the subset.
#' @return The mean pairwise patristic distance.
mean_pd <- function(dist_mat, subset) {
  if (length(subset) < 2) return(0)
  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  upper_vals <- sub_dist[upper.tri(sub_dist)]
  mean(upper_vals, na.rm = TRUE)
}

#' @rdname distances
#' @export
#' @title Mean nearest-neighbor distance
#' @param dist_mat A symmetric patristic distance matrix.
#' @param subset Character vector of species names in the subset.
#' @return The mean nearest-neighbor phylogenetic distance.
mean_nnd <- function(dist_mat, subset) {
  if (length(subset) < 2) return(0)
  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  diag(sub_dist) <- Inf
  nnd <- apply(sub_dist, 1, min, na.rm = TRUE)
  mean(nnd, na.rm = TRUE)
}

#' @rdname distances
#' @export
#' @title Nearest-neighbor distances
#' @description Returns the nearest-neighbor distance for each species in the
#'   subset.
#' @param dist_mat A symmetric patristic distance matrix.
#' @param subset Character vector of species names in the subset.
#' @return A named numeric vector of nearest-neighbor distances.
nearest_neighbor_distances <- function(dist_mat, subset) {
  if (length(subset) < 2) {
    return(stats::setNames(0, subset))
  }
  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  diag(sub_dist) <- Inf
  nnd <- apply(sub_dist, 1, min, na.rm = TRUE)
  nnd
}
