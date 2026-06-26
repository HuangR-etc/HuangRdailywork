#' Distance-based subset descriptors
#'
#' Implements the within-subset phylogenetic distance criteria used in the
#' manuscript: MinPD, MeanPD, MeanNND, and MaxPD.
NULL

distance_metrics <- function(dist_mat, subset) {
  if (length(subset) < 2) {
    return(list(MinPD = 0, MeanPD = 0, MeanNND = 0, MaxPD = 0))
  }

  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  diag(sub_dist) <- NA_real_
  upper_vals <- sub_dist[upper.tri(sub_dist)]

  sub_dist_nnd <- sub_dist
  diag(sub_dist_nnd) <- Inf
  nnd <- apply(sub_dist_nnd, 1, min, na.rm = TRUE)

  list(
    MinPD = min(upper_vals, na.rm = TRUE),
    MeanPD = mean(upper_vals, na.rm = TRUE),
    MeanNND = mean(nnd, na.rm = TRUE),
    MaxPD = max(upper_vals, na.rm = TRUE)
  )
}

#' @export
min_pd <- function(dist_mat, subset) distance_metrics(dist_mat, subset)$MinPD

#' @export
mean_pd <- function(dist_mat, subset) distance_metrics(dist_mat, subset)$MeanPD

#' @export
mean_nnd <- function(dist_mat, subset) distance_metrics(dist_mat, subset)$MeanNND

#' @export
max_pd <- function(dist_mat, subset) distance_metrics(dist_mat, subset)$MaxPD

#' @export
nearest_neighbor_distances <- function(dist_mat, subset) {
  if (length(subset) < 2) {
    return(stats::setNames(rep(0, length(subset)), subset))
  }

  sub_dist <- dist_mat[subset, subset, drop = FALSE]
  diag(sub_dist) <- Inf
  apply(sub_dist, 1, min, na.rm = TRUE)
}

subset_metrics_table <- function(dist_mat, subsets) {
  do.call(rbind, lapply(seq_along(subsets), function(i) {
    m <- distance_metrics(dist_mat, subsets[[i]])
    data.frame(
      SubsetID = i,
      MinPD = m$MinPD,
      MeanPD = m$MeanPD,
      MeanNND = m$MeanNND,
      MaxPD = m$MaxPD,
      stringsAsFactors = FALSE
    )
  }))
}
