#' Dispersed subset selection
#'
#' Selects a phylogenetically dispersed subset using greedy forward selection
#' followed by optional swap-based local refinement. The objective prioritizes
#' maximizing MinPD, then MeanPD, then MeanNND (lexicographic order).
#'
#' @name select_dispersed
NULL

#' Lexicographic comparison for maximization
#' @param a List with MinPD, MeanPD, MeanNND
#' @param b List with MinPD, MeanPD, MeanNND
#' @return TRUE if a is better than b
#' @keywords internal
is_better_lexico_max <- function(a, b) {
  if (a$MinPD > b$MinPD) return(TRUE)
  if (a$MinPD < b$MinPD) return(FALSE)
  if (a$MeanPD > b$MeanPD) return(TRUE)
  if (a$MeanPD < b$MeanPD) return(FALSE)
  if (a$MeanNND > b$MeanNND) return(TRUE)
  FALSE
}

#' Find the species with maximum mean distance to all others
#' @param dist_mat Distance matrix
#' @return List with index, name, and mean distance
#' @keywords internal
find_peripheral_species <- function(dist_mat) {
  d <- dist_mat
  diag(d) <- NA
  mean_dist <- rowMeans(d, na.rm = TRUE)
  idx <- which.max(mean_dist)
  list(idx = idx, name = names(idx), mean_distance = mean_dist[idx])
}

#' Greedy forward selection for dispersed subset
#' @param dist_mat Distance matrix
#' @param candidates Character vector of candidate species
#' @param size Target subset size
#' @param start Starting strategy: "peripheral" or "random"
#' @param start_subset Optional user-provided starting subset
#' @param seed Random seed for "random" start
#' @return List with selected, greedy_path, metrics
#' @keywords internal
greedy_dispersed <- function(dist_mat, candidates, size,
                             start = c("peripheral", "random", "user"),
                             start_subset = NULL, seed = NULL) {
  start <- match.arg(start)

  if (start == "user") {
    if (is.null(start_subset)) {
      stop("start_subset must be provided when start = 'user'")
    }
    current <- start_subset
  } else if (start == "random") {
    if (!is.null(seed)) set.seed(seed)
    current <- sample(candidates, 1)
  } else {
    # peripheral: species with highest mean distance
    peripheral <- find_peripheral_species(dist_mat)
    current <- peripheral$name
  }

  available <- setdiff(candidates, current)
  greedy_path <- list(current)

  while (length(current) < size) {
    best_candidate <- NULL
    best_metrics <- NULL

    for (cand in available) {
      temp <- c(current, cand)
      temp_metrics <- distance_metrics(dist_mat, temp)

      if (is.null(best_candidate) ||
          is_better_lexico_max(temp_metrics, best_metrics)) {
        best_candidate <- cand
        best_metrics <- temp_metrics
      }
    }

    current <- c(current, best_candidate)
    available <- setdiff(available, best_candidate)
    greedy_path[[length(greedy_path) + 1]] <- current
  }

  list(
    selected = current,
    greedy_path = greedy_path,
    metrics = distance_metrics(dist_mat, current)
  )
}

#' Swap-based local refinement for dispersed subset
#' @param dist_mat Distance matrix
#' @param current Character vector of current subset
#' @param candidates Full candidate pool
#' @param max_iter Maximum iterations
#' @return List with selected, metrics, n_swaps, converged
#' @keywords internal
swap_refine <- function(dist_mat, current, candidates,
                        max_iter = Inf) {
  n <- length(current)
  available <- setdiff(candidates, current)
  current_metrics <- distance_metrics(dist_mat, current)
  improved <- TRUE
  n_iter <- 0
  n_swaps <- 0

  while (improved && n_iter < max_iter) {
    improved <- FALSE
    n_iter <- n_iter + 1

    for (out_sp in current) {
      for (in_sp in available) {
        new_subset <- setdiff(current, out_sp)
        new_subset <- c(new_subset, in_sp)
        new_metrics <- distance_metrics(dist_mat, new_subset)

        if (is_better_lexico_max(new_metrics, current_metrics)) {
          current <- new_subset
          current_metrics <- new_metrics
          available <- c(available, out_sp)
          available <- setdiff(available, in_sp)
          improved <- TRUE
          n_swaps <- n_swaps + 1
          break
        }
      }
      if (improved) break
    }
  }

  list(
    selected = current,
    metrics = current_metrics,
    n_swaps = n_swaps,
    n_iter = n_iter,
    converged = !improved
  )
}

#' @rdname select_dispersed
#' @export
#' @title Select a phylogenetically dispersed subset
#' @description Selects a dispersed subset from candidate species using greedy
#'   forward selection followed by optional swap-based local refinement.
#'   The objective is lexicographic: maximize MinPD > MeanPD > MeanNND.
#' @param tree Optional phylogenetic tree of class \code{phylo}. If provided,
#'   the patristic distance matrix is computed automatically.
#' @param dist_mat Optional patristic distance matrix. If \code{tree} is
#'   provided, this is ignored.
#' @param candidates Character vector of candidate species names.
#' @param size Target subset size.
#' @param objective Objective ordering. Currently only \code{"MinPD"} is
#'   supported (lexicographic: MinPD > MeanPD > MeanNND).
#' @param start Starting strategy: \code{"peripheral"} (species with highest
#'   mean distance), \code{"random"}, or \code{"user"}.
#' @param start_subset Optional starting subset when \code{start = "user"}.
#' @param refine If \code{TRUE}, perform swap-based local refinement after
#'   greedy selection.
#' @param max_iter Maximum iterations for local refinement. Default \code{Inf}
#'   means run until convergence.
#' @param seed Random seed for reproducibility.
#' @return A list with components:
#'   \item{selected}{Character vector of selected species.}
#'   \item{greedy_path}{List of subsets at each greedy step.}
#'   \item{swaps}{List of swaps accepted during refinement (if refine=TRUE).}
#'   \item{metrics}{Distance metrics for the final subset.}
#'   \item{converged}{Logical, whether local refinement reached a local optimum.}
#'   \item{call}{The matched function call.}
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(20)
#' candidates <- tree$tip.label
#' result <- select_dispersed(tree = tree, candidates = candidates, size = 5)
#' result$selected
#' }
select_dispersed <- function(tree = NULL, dist_mat = NULL, candidates, size,
                             objective = c("MinPD", "MeanPD", "MeanNND"),
                             start = c("peripheral", "random", "user"),
                             start_subset = NULL,
                             refine = TRUE, max_iter = Inf, seed = NULL) {
  objective <- match.arg(objective)
  start <- match.arg(start)

  if (is.null(dist_mat)) {
    if (is.null(tree)) {
      stop("Either tree or dist_mat must be provided")
    }
    dist_mat <- patristic_matrix(tree)
  }

  if (!is.null(seed)) set.seed(seed)

  # Greedy selection
  greedy_result <- greedy_dispersed(
    dist_mat = dist_mat,
    candidates = candidates,
    size = size,
    start = start,
    start_subset = start_subset,
    seed = seed
  )

  swaps <- NULL
  converged <- NULL

  if (refine && size >= 2) {
    refine_result <- swap_refine(
      dist_mat = dist_mat,
      current = greedy_result$selected,
      candidates = candidates,
      max_iter = max_iter
    )
    final_selected <- refine_result$selected
    final_metrics <- refine_result$metrics
    converged <- refine_result$converged
    swaps <- refine_result$n_swaps
  } else {
    final_selected <- greedy_result$selected
    final_metrics <- greedy_result$metrics
    converged <- NA
  }

  structure(list(
    selected = final_selected,
    greedy_path = greedy_result$greedy_path,
    swaps = swaps,
    metrics = final_metrics,
    converged = converged,
    call = match.call()
  ), class = "phylosubset_selection")
}
