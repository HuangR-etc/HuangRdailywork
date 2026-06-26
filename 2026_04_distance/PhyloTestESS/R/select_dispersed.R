#' Dispersed subset selection
#'
#' The manuscript default is a two-stage heuristic:
#' 1. greedy forward selection starting from the most peripheral species
#' 2. one-for-one exchange refinement under the same lexicographic objective
#'
#' Objective: maximize MinPD, then MeanPD, then MeanNND.
NULL

is_better_dispersed <- function(a, b, tol = 1e-10) {
  if (is.null(b)) return(TRUE)
  if (a$MinPD > b$MinPD + tol) return(TRUE)
  if (a$MinPD < b$MinPD - tol) return(FALSE)
  if (a$MeanPD > b$MeanPD + tol) return(TRUE)
  if (a$MeanPD < b$MeanPD - tol) return(FALSE)
  if (a$MeanNND > b$MeanNND + tol) return(TRUE)
  FALSE
}

peripheral_species <- function(dist_mat, candidates) {
  d <- dist_mat[candidates, candidates, drop = FALSE]
  diag(d) <- NA_real_
  mean_dist <- rowMeans(d, na.rm = TRUE)
  names(which.max(mean_dist))
}

greedy_dispersed <- function(dist_mat, candidates, size, start_species = NULL, tol = 1e-10) {
  if (is.null(start_species)) {
    start_species <- peripheral_species(dist_mat, candidates)
  }

  selected <- start_species
  greedy_path <- list(selected)

  while (length(selected) < size) {
    available <- setdiff(candidates, selected)
    best_species <- NULL
    best_metrics <- NULL

    for (cand in available) {
      metrics <- distance_metrics(dist_mat, c(selected, cand))
      if (is_better_dispersed(metrics, best_metrics, tol = tol)) {
        best_species <- cand
        best_metrics <- metrics
      }
    }

    selected <- c(selected, best_species)
    greedy_path[[length(greedy_path) + 1]] <- selected
  }

  list(
    selected = selected,
    greedy_path = greedy_path,
    metrics = distance_metrics(dist_mat, selected),
    start_species = start_species
  )
}

swap_refine_dispersed <- function(dist_mat, selected, candidates, max_iter = Inf, tol = 1e-10) {
  current <- selected
  current_metrics <- distance_metrics(dist_mat, current)
  swap_log <- list()
  iteration <- 0
  converged <- FALSE

  repeat {
    if (iteration >= max_iter) break

    improved <- FALSE
    available <- setdiff(candidates, current)

    for (out_sp in current) {
      for (in_sp in available) {
        proposal <- c(setdiff(current, out_sp), in_sp)
        proposal_metrics <- distance_metrics(dist_mat, proposal)

        if (is_better_dispersed(proposal_metrics, current_metrics, tol = tol)) {
          iteration <- iteration + 1
          current <- proposal
          current_metrics <- proposal_metrics
          swap_log[[iteration]] <- list(
            iteration = iteration,
            out_species = out_sp,
            in_species = in_sp,
            metrics = proposal_metrics
          )
          improved <- TRUE
          break
        }
      }
      if (improved) break
    }

    if (!improved) {
      converged <- TRUE
      break
    }
  }

  list(
    selected = current,
    metrics = current_metrics,
    swap_log = swap_log,
    iterations = iteration,
    converged = converged
  )
}

#' @export
select_dispersed <- function(tree = NULL, dist_mat = NULL, candidates, size,
                             refine = TRUE, max_iter = Inf, tol = 1e-10) {
  if (is.null(dist_mat)) {
    if (is.null(tree)) {
      stop("Either tree or dist_mat must be provided.", call. = FALSE)
    }
    dist_mat <- patristic_matrix(tree, candidates)
  }

  if (size >= length(candidates)) {
    stop("size must be smaller than the number of candidates.", call. = FALSE)
  }

  greedy <- greedy_dispersed(dist_mat, candidates, size, tol = tol)

  if (!refine || size < 2) {
    return(list(
      selected = greedy$selected,
      metrics = greedy$metrics,
      start_species = greedy$start_species,
      greedy_path = greedy$greedy_path,
      swap_log = list(),
      converged = NA,
      algorithm = "dispersed_fixed_max_mean_distance_greedy"
    ))
  }

  refined <- swap_refine_dispersed(
    dist_mat = dist_mat,
    selected = greedy$selected,
    candidates = candidates,
    max_iter = max_iter,
    tol = tol
  )

  list(
    selected = refined$selected,
    metrics = refined$metrics,
    start_species = greedy$start_species,
    greedy_path = greedy$greedy_path,
    swap_log = refined$swap_log,
    iterations = refined$iterations,
    converged = refined$converged,
    algorithm = "dispersed_fixed_start_greedy_exchange"
  )
}
