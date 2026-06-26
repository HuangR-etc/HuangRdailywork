#' Main package entry points
NULL

resolve_subset_selection <- function(tree, dist_mat, candidates, size,
                                     subset_type = c("dispersed", "clustered"),
                                     clustered_method = c("multistart_exchange", "fast")) {
  subset_type <- match.arg(subset_type)
  clustered_method <- match.arg(clustered_method)

  if (subset_type == "dispersed") {
    return(select_dispersed(tree = tree, dist_mat = dist_mat, candidates = candidates, size = size))
  }

  select_clustered(
    tree = tree,
    dist_mat = dist_mat,
    candidates = candidates,
    size = size,
    method = clustered_method
  )
}

#' @export
phylo_test_ess <- function(tree, candidates, size,
                           subset_type = c("dispersed", "clustered"),
                           clustered_method = c("multistart_exchange", "fast"),
                           n_random = 1000,
                           cov_model = c("BM", "lambda", "OU", "EB"),
                           lambda = 1,
                           half_life = NULL,
                           alpha = NULL,
                           eb_rate = NULL,
                           eb_method = c("simple", "branch_transform"),
                           compute_piess = FALSE,
                           piess_n_sim = 10000,
                           piess_error_sd = sqrt(0.1),
                           seed = NULL) {
  subset_type <- match.arg(subset_type)
  clustered_method <- match.arg(clustered_method)
  cov_model <- match.arg(cov_model)
  eb_method <- match.arg(eb_method)

  check_phylo_input(tree, candidates)
  if (size >= length(candidates)) {
    stop("size must be smaller than the number of candidates.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  dist_mat <- patristic_matrix(tree, candidates)
  selection <- resolve_subset_selection(
    tree = tree,
    dist_mat = dist_mat,
    candidates = candidates,
    size = size,
    subset_type = subset_type,
    clustered_method = clustered_method
  )

  observed_metrics <- selection$metrics
  random_metrics <- random_baseline(
    dist_mat = dist_mat,
    candidates = candidates,
    size = size,
    n = n_random,
    seed = if (!is.null(seed)) seed + 1 else NULL
  )
  random_comparison <- compare_to_random(observed_metrics, random_metrics, type = subset_type)

  V <- phylo_covariance(
    tree = tree,
    tips = candidates,
    model = cov_model,
    lambda = lambda,
    half_life = half_life,
    alpha = alpha,
    eb_rate = eb_rate,
    eb_method = eb_method
  )

  R_sub <- cov_to_cor(V[selection$selected, selection$selected, drop = FALSE])
  dep <- dependence_metrics(R_sub)

  piess <- NULL
  if (compute_piess) {
    piess <- analyze_prediction_metric_ess(
      R_full = cov_to_cor(V),
      subset_names = selection$selected,
      subset_type = subset_type,
      n_sim = piess_n_sim,
      error_sd = piess_error_sd,
      seed = if (!is.null(seed)) seed + 2 else PRED_ESS_SEED
    )
  }

  structure(
    list(
      call = match.call(),
      package = "PhyloTestESS",
      subset_type = subset_type,
      clustered_method = if (subset_type == "clustered") clustered_method else NA_character_,
      selected = selection$selected,
      size = size,
      candidates = candidates,
      selection = selection,
      distance_metrics = observed_metrics,
      random_metrics = random_metrics,
      random_comparison = random_comparison,
      covariance_model = cov_model,
      covariance_matrix = V,
      correlation_matrix = R_sub,
      dependence_metrics = dep,
      piess = piess
    ),
    class = "phylo_test_ess_result"
  )
}

#' @export
phylo_subset <- function(tree, candidates, size,
                         type = c("dispersed", "clustered"),
                         clustered_method = c("multistart_exchange", "fast"),
                         n_random = 1000,
                         cov_model = c("BM", "lambda", "OU", "EB"),
                         lambda = 1,
                         half_life = NULL,
                         alpha = NULL,
                         eb_rate = NULL,
                         eb_method = c("simple", "branch_transform"),
                         compute_piess = FALSE,
                         piess_n_sim = 10000,
                         piess_error_sd = sqrt(0.1),
                         seed = NULL) {
  type <- match.arg(type)
  clustered_method <- match.arg(clustered_method)
  cov_model <- match.arg(cov_model)
  eb_method <- match.arg(eb_method)

  phylo_test_ess(
    tree = tree,
    candidates = candidates,
    size = size,
    subset_type = type,
    clustered_method = clustered_method,
    n_random = n_random,
    cov_model = cov_model,
    lambda = lambda,
    half_life = half_life,
    alpha = alpha,
    eb_rate = eb_rate,
    eb_method = eb_method,
    compute_piess = compute_piess,
    piess_n_sim = piess_n_sim,
    piess_error_sd = piess_error_sd,
    seed = seed
  )
}

#' @export
print.phylo_test_ess_result <- function(x, ...) {
  cat("PhyloTestESS result\n")
  cat("Subset type:", x$subset_type, "\n")
  cat("Selected species:", length(x$selected), "of", length(x$candidates), "\n")
  if (!is.na(x$clustered_method)) {
    cat("Clustered method:", x$clustered_method, "\n")
  }
  cat("Distance metrics:\n")
  print(round(as.data.frame(x$distance_metrics), 4))
  cat("Dependence metrics:\n")
  print(round(as.data.frame(x$dependence_metrics), 4))
  invisible(x)
}

#' @export
summary.phylo_test_ess_result <- function(object, ...) {
  print(object)
  cat("Random comparison:\n")
  print(object$random_comparison)
  if (!is.null(object$piess)) {
    cat("PIESS:\n")
    print(object$piess$piess)
  }
  invisible(object)
}

#' @export
as.data.frame.phylo_test_ess_result <- function(x, ...) {
  data.frame(
    subset_type = x$subset_type,
    size = x$size,
    n_candidates = length(x$candidates),
    n_selected = length(x$selected),
    MinPD = x$distance_metrics$MinPD,
    MeanPD = x$distance_metrics$MeanPD,
    MeanNND = x$distance_metrics$MeanNND,
    MaxPD = x$distance_metrics$MaxPD,
    MeanOffCor = x$dependence_metrics$MeanOffCor,
    MaxOffCor = x$dependence_metrics$MaxOffCor,
    MIESS = x$dependence_metrics$MIESS,
    stringsAsFactors = FALSE
  )
}
