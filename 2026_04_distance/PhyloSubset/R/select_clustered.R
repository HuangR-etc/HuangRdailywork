#' Clustered (high-dependence) subset selection
#'
#' Constructs a clustered reference subset by taking each species as a seed,
#' finding its nearest (size - 1) neighbors, and selecting the neighborhood
#' with the smallest MeanPD, MeanNND, and MaxPD.
#'
#' @name select_clustered
NULL

#' Generate clustered candidate subsets via seed-neighborhood
#' @param dist_mat Distance matrix
#' @param candidates Character vector of candidate species
#' @param size Target subset size
#' @return List with subsets, seed_idx, seed_names
#' @keywords internal
generate_seed_neighborhoods <- function(dist_mat, candidates, size) {
  d <- dist_mat[candidates, candidates, drop = FALSE]

  subsets <- vector("list", length(candidates))
  names(subsets) <- candidates

  for (seed in candidates) {
    distances_from_seed <- d[seed, ]
    distances_from_seed[seed] <- Inf

    nearest <- names(sort(distances_from_seed)[1:(size - 1)])
    subset_now <- sort(c(seed, nearest))
    subsets[[seed]] <- subset_now
  }

  # Remove duplicates
  keys <- vapply(subsets, function(x) paste(x, collapse = "|"), character(1))
  unique_keys <- unique(keys)
  unique_subsets <- subsets[match(unique_keys, keys)]
  unique_seeds <- names(unique_subsets)

  list(
    subsets = unname(unique_subsets),
    seed_names = unique_seeds,
    n_raw = length(subsets),
    n_unique = length(unique_subsets)
  )
}

#' @rdname select_clustered
#' @export
#' @title Select a phylogenetically clustered (high-dependence) subset
#' @description Constructs a clustered reference subset by taking each species
#'   as a seed, finding its nearest (size - 1) neighbors, and selecting the
#'   neighborhood with the smallest MeanPD, MeanNND, and MaxPD. This is
#'   intended as a high-dependence comparison group, not a primary optimization
#'   target.
#' @param tree Optional phylogenetic tree of class \code{phylo}. If provided,
#'   the patristic distance matrix is computed automatically.
#' @param dist_mat Optional patristic distance matrix.
#' @param candidates Character vector of candidate species names.
#' @param size Target subset size.
#' @param objective Objective ordering: \code{"MeanPD"} (default), \code{"MeanNND"},
#'   or \code{"MaxPD"}.
#' @param seed_species Optional character vector of seed species to use. If
#'   NULL, all candidates are used as seeds.
#' @param collapse_duplicates If \code{TRUE} (default), remove duplicate
#'   neighborhoods (same set of species from different seeds).
#' @return A list with components:
#'   \item{selected}{Character vector of selected species.}
#'   \item{seed}{The seed species that generated the best neighborhood.}
#'   \item{metrics}{Distance metrics for the selected subset.}
#'   \item{candidate_metrics}{Data frame of metrics for all candidate
#'     neighborhoods.}
#'   \item{n_candidates}{Number of unique neighborhoods evaluated.}
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(20)
#' candidates <- tree$tip.label
#' result <- select_clustered(tree = tree, candidates = candidates, size = 5)
#' result$selected
#' }
select_clustered <- function(tree = NULL, dist_mat = NULL, candidates, size,
                             objective = c("MeanPD", "MeanNND", "MaxPD"),
                             seed_species = NULL,
                             collapse_duplicates = TRUE) {
  objective <- match.arg(objective)

  if (is.null(dist_mat)) {
    if (is.null(tree)) {
      stop("Either tree or dist_mat must be provided")
    }
    dist_mat <- patristic_matrix(tree)
  }

  if (!is.null(seed_species)) {
    candidates <- intersect(seed_species, candidates)
    if (length(candidates) == 0) {
      stop("No seed species found in candidates")
    }
  }

  # Generate neighborhoods
  neighborhoods <- generate_seed_neighborhoods(
    dist_mat = dist_mat,
    candidates = candidates,
    size = size
  )

  # Compute metrics for each neighborhood
  metrics_list <- lapply(neighborhoods$subsets, function(sub) {
    distance_metrics(dist_mat, sub)
  })

  metrics_df <- do.call(rbind, lapply(metrics_list, function(m) {
    data.frame(
      MinPD = m$MinPD,
      MeanPD = m$MeanPD,
      MeanNND = m$MeanNND,
      MaxPD = m$MaxPD,
      stringsAsFactors = FALSE
    )
  }))

  metrics_df$Seed <- neighborhoods$seed_names

  # Sort by objective (ascending)
  ord <- order(metrics_df[[objective]], metrics_df$MeanNND, metrics_df$MaxPD)
  best_idx <- ord[1]

  list(
    selected = neighborhoods$subsets[[best_idx]],
    seed = neighborhoods$seed_names[best_idx],
    metrics = metrics_list[[best_idx]],
    candidate_metrics = metrics_df[ord, , drop = FALSE],
    n_candidates = neighborhoods$n_unique
  )
}
