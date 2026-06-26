#' Clustered subset selection
#'
#' Main method:
#' multistart greedy forward selection minimizing MeanPD, MeanNND, MaxPD,
#' followed by one-for-one exchange refinement.
#'
#' Fast fallback:
#' seed-based nearest-neighbor expansion retained from the old package.
NULL

is_better_clustered <- function(a, b, tol = 1e-10) {
  if (is.null(b)) return(TRUE)
  if (a$MeanPD < b$MeanPD - tol) return(TRUE)
  if (a$MeanPD > b$MeanPD + tol) return(FALSE)
  if (a$MeanNND < b$MeanNND - tol) return(TRUE)
  if (a$MeanNND > b$MeanNND + tol) return(FALSE)
  if (a$MaxPD < b$MaxPD - tol) return(TRUE)
  FALSE
}

greedy_clustered_from_start <- function(dist_mat, candidates, size, start_species, tol = 1e-10) {
  selected <- start_species
  greedy_log <- list()

  while (length(selected) < size) {
    available <- setdiff(candidates, selected)
    best_species <- NULL
    best_metrics <- NULL

    for (cand in available) {
      metrics <- distance_metrics(dist_mat, c(selected, cand))
      if (is_better_clustered(metrics, best_metrics, tol = tol)) {
        best_species <- cand
        best_metrics <- metrics
      }
    }

    selected <- c(selected, best_species)
    greedy_log[[length(greedy_log) + 1]] <- list(
      step = length(selected),
      added_species = best_species,
      metrics = best_metrics
    )
  }

  list(
    selected = selected,
    metrics = distance_metrics(dist_mat, selected),
    greedy_log = greedy_log,
    start_species = start_species
  )
}

refine_clustered_exchange <- function(dist_mat, selected, candidates,
                                      max_iter = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
                                      tol = CLUSTERED_EXCHANGE_TOL) {
  current <- selected
  current_metrics <- distance_metrics(dist_mat, current)
  exchange_log <- list()
  iteration <- 0
  converged <- FALSE

  repeat {
    if (iteration >= max_iter) break
    available <- setdiff(candidates, current)
    best_swap <- NULL
    best_subset <- NULL
    best_metrics <- NULL

    for (out_sp in current) {
      for (in_sp in available) {
        proposal <- c(setdiff(current, out_sp), in_sp)
        proposal_metrics <- distance_metrics(dist_mat, proposal)

        if (is_better_clustered(proposal_metrics, current_metrics, tol = tol) &&
            is_better_clustered(proposal_metrics, best_metrics, tol = tol)) {
          best_swap <- c(out_sp, in_sp)
          best_subset <- proposal
          best_metrics <- proposal_metrics
        }
      }
    }

    if (is.null(best_subset)) {
      converged <- TRUE
      break
    }

    iteration <- iteration + 1
    current <- best_subset
    current_metrics <- best_metrics
    exchange_log[[iteration]] <- list(
      iteration = iteration,
      out_species = best_swap[[1]],
      in_species = best_swap[[2]],
      metrics = best_metrics
    )
  }

  list(
    selected = current,
    metrics = current_metrics,
    exchange_log = exchange_log,
    iterations = iteration,
    converged = converged
  )
}

clustered_seed_neighborhoods <- function(dist_mat, candidates, size) {
  subsets <- lapply(candidates, function(seed) {
    distances <- dist_mat[seed, candidates]
    distances[seed] <- Inf
    neighbors <- names(sort(distances))[seq_len(size - 1)]
    sort(c(seed, neighbors))
  })

  keys <- vapply(subsets, paste, collapse = "|", character(1))
  unique_keys <- unique(keys)
  unique_subsets <- subsets[match(unique_keys, keys)]
  unique_seeds <- candidates[match(unique_keys, keys)]

  list(
    subsets = unique_subsets,
    seeds = unique_seeds,
    n_raw = length(subsets),
    n_unique = length(unique_subsets)
  )
}

#' @export
select_clustered_fast <- function(tree = NULL, dist_mat = NULL, candidates, size) {
  if (is.null(dist_mat)) {
    if (is.null(tree)) {
      stop("Either tree or dist_mat must be provided.", call. = FALSE)
    }
    dist_mat <- patristic_matrix(tree, candidates)
  }

  cand <- clustered_seed_neighborhoods(dist_mat, candidates, size)
  metrics <- subset_metrics_table(dist_mat, cand$subsets)
  metrics$Seed_Name <- cand$seeds

  ord <- order(metrics$MeanPD, metrics$MeanNND, metrics$MaxPD)
  best_idx <- ord[[1]]

  list(
    selected = cand$subsets[[best_idx]],
    metrics = distance_metrics(dist_mat, cand$subsets[[best_idx]]),
    candidate_metrics = metrics,
    ordered_candidate_metrics = metrics[ord, , drop = FALSE],
    seed = cand$seeds[[best_idx]],
    best_seed_name = cand$seeds[[best_idx]],
    n_raw_candidates = cand$n_raw,
    n_unique_candidates = cand$n_unique,
    algorithm = "clustered_seed_nearest_neighbors_meanpd_meannnd_maxpd"
  )
}

#' @export
select_clustered <- function(tree = NULL, dist_mat = NULL, candidates, size,
                             method = c("multistart_exchange", "fast"),
                             max_iter = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
                             tol = CLUSTERED_EXCHANGE_TOL,
                             verbose = FALSE) {
  method <- match.arg(method)

  if (is.null(dist_mat)) {
    if (is.null(tree)) {
      stop("Either tree or dist_mat must be provided.", call. = FALSE)
    }
    dist_mat <- patristic_matrix(tree, candidates)
  }

  if (size >= length(candidates)) {
    stop("size must be smaller than the number of candidates.", call. = FALSE)
  }

  if (method == "fast") {
    return(select_clustered_fast(
      dist_mat = dist_mat,
      candidates = candidates,
      size = size
    ))
  }

  start_results <- vector("list", length(candidates))

  for (i in seq_along(candidates)) {
    start_species <- candidates[[i]]
    if (verbose) {
      message("Evaluating clustered start ", i, "/", length(candidates), ": ", start_species)
    }

    greedy <- greedy_clustered_from_start(
      dist_mat = dist_mat,
      candidates = candidates,
      size = size,
      start_species = start_species,
      tol = tol
    )

    refined <- refine_clustered_exchange(
      dist_mat = dist_mat,
      selected = greedy$selected,
      candidates = candidates,
      max_iter = max_iter,
      tol = tol
    )

    start_results[[i]] <- list(
      start_species = start_species,
      greedy = greedy,
      refined = refined,
      selected = refined$selected,
      metrics = refined$metrics
    )
  }

  ranking <- do.call(rbind, lapply(seq_along(start_results), function(i) {
    m <- start_results[[i]]$metrics
    gm <- start_results[[i]]$greedy$metrics
    data.frame(
      Start_Name = start_results[[i]]$start_species,
      Final_MinPD = m$MinPD,
      Final_MeanPD = m$MeanPD,
      Final_MeanNND = m$MeanNND,
      Final_MaxPD = m$MaxPD,
      Greedy_MinPD = gm$MinPD,
      Greedy_MeanPD = gm$MeanPD,
      Greedy_MeanNND = gm$MeanNND,
      Greedy_MaxPD = gm$MaxPD,
      Exchange_Converged = start_results[[i]]$refined$converged,
      MinPD = m$MinPD,
      MeanPD = m$MeanPD,
      MeanNND = m$MeanNND,
      MaxPD = m$MaxPD,
      Exchange_Iterations = start_results[[i]]$refined$iterations,
      stringsAsFactors = FALSE
    )
  }))

  ranking$Final_Subset_Key <- vapply(
    start_results,
    function(x) paste(sort(x$selected), collapse = "|"),
    character(1)
  )

  ord <- order(ranking$MeanPD, ranking$MeanNND, ranking$MaxPD)
  best <- start_results[[ord[[1]]]]

  list(
    selected = best$selected,
    metrics = best$metrics,
    candidate_metrics = ranking,
    ordered_candidate_metrics = ranking[ord, , drop = FALSE],
    start_results = start_results,
    best_start = best$start_species,
    best_start_name = best$start_species,
    n_raw_candidates = length(candidates),
    n_unique_candidates = length(unique(ranking$Final_Subset_Key)),
    algorithm = "clustered_multistart_greedy_exchange_meanpd_meannnd_maxpd"
  )
}
