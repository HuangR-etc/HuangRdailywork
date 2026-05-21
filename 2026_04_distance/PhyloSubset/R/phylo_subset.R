#' PhyloSubset: Main entry point and S3 class
#'
#' The main user-facing function that runs the complete analysis pipeline:
#' subset selection, random baseline, empirical p-values, and dependence
#' diagnostics.
#'
#' @name phylo_subset
NULL

#' @rdname phylo_subset
#' @export
#' @title One-step phylogenetic subset analysis
#' @description The main entry point for PhyloSubset. Given a phylogenetic
#'   tree, a candidate species pool, and a target subset size, this function
#'   automatically selects a dispersed or clustered subset, computes distance
#'   metrics (MinPD, MeanPD, MeanNND, MaxPD), generates a random baseline,
#'   calculates empirical p-values, and evaluates phylogenetic dependence
#'   diagnostics (MeanOffCor, MaxOffCor, MeanESS).
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param candidates Character vector of candidate species names. Must be
#'   present in \code{tree}.
#' @param size Target subset size.
#' @param type Subset type: \code{"dispersed"} (default) or \code{"clustered"}.
#' @param n_random Number of random subsets for baseline (default: 1000).
#' @param cov_model Covariance model for dependence diagnostics:
#'   \code{"BM"} (default), \code{"lambda"}, or \code{"OU"}.
#' @param lambda Lambda parameter for Pagel's lambda model (default: 1).
#'   Ignored unless \code{cov_model = "lambda"}.
#' @param half_life Half-life parameter for OU model. Ignored unless
#'   \code{cov_model = "OU"}.
#' @param seed Optional random seed for reproducibility.
#' @param return_random If \code{TRUE} (default), include random baseline
#'   metrics in the result.
#' @return An object of class \code{phylosubset_result} containing:
#'   \item{call}{The matched function call.}
#'   \item{type}{Subset type.}
#'   \item{selected}{Character vector of selected species.}
#'   \item{candidates}{Character vector of candidate species.}
#'   \item{size}{Target subset size.}
#'   \item{distance_metrics}{Data frame of observed distance metrics.}
#'   \item{dependence_metrics}{Data frame of dependence diagnostics.}
#'   \item{random_metrics}{Data frame of random baseline metrics (if
#'     \code{return_random = TRUE}).}
#'   \item{empirical_p}{Data frame of empirical p-values and SES.}
#'   \item{tree}{The input tree.}
#'   \item{dist_mat}{The patristic distance matrix.}
#'   \item{settings}{List of analysis settings.}
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(30)
#' candidates <- tree$tip.label
#'
#' # Basic usage
#' res <- phylo_subset(tree, candidates, size = 8, n_random = 100, seed = 123)
#' summary(res)
#'
#' # Clustered subset
#' res_clust <- phylo_subset(tree, candidates, size = 8, type = "clustered",
#'                           n_random = 100, seed = 123)
#' summary(res_clust)
#' }
phylo_subset <- function(tree, candidates, size,
                         type = c("dispersed", "clustered"),
                         n_random = 1000,
                         cov_model = c("BM", "lambda", "OU"),
                         lambda = 1,
                         half_life = NULL,
                         seed = NULL,
                         return_random = TRUE) {
  type <- match.arg(type)
  cov_model <- match.arg(cov_model)

  # Validate inputs
  check_phylo_input(tree, candidates)

  if (size >= length(candidates)) {
    stop("size must be smaller than the number of candidates")
  }

  if (!is.null(seed)) set.seed(seed)

  # Compute distance matrix
  dist_mat <- patristic_matrix(tree, candidates)

  # Select subset
  if (type == "dispersed") {
    sel <- select_dispersed(
      dist_mat = dist_mat,
      candidates = candidates,
      size = size,
      refine = TRUE,
      seed = seed
    )
  } else {
    sel <- select_clustered(
      dist_mat = dist_mat,
      candidates = candidates,
      size = size
    )
  }

  # Observed distance metrics
  obs_metrics <- distance_metrics(dist_mat, sel$selected)

  # Random baseline
  if (return_random || TRUE) {
    rand_metrics <- random_baseline(
      dist_mat = dist_mat,
      candidates = candidates,
      size = size,
      n = n_random,
      seed = if (!is.null(seed)) seed + 1 else NULL
    )
  } else {
    rand_metrics <- NULL
  }

  # Empirical p-values
  emp_p <- compare_to_random(obs_metrics, rand_metrics, type = type)

  # Covariance and dependence
  V <- phylo_covariance(
    tree = tree,
    tips = candidates,
    model = cov_model,
    lambda = lambda,
    half_life = half_life
  )

  R_sub <- cov_to_cor(V[sel$selected, sel$selected, drop = FALSE])
  dep <- dependence_metrics(R_sub)

  # Random dependence
  if (return_random) {
    rand_subsets <- random_subsets(candidates, size, n_random,
                                   seed = if (!is.null(seed)) seed + 2 else NULL)
    rand_dep <- do.call(rbind, lapply(rand_subsets, function(sub) {
      R_rand <- cov_to_cor(V[sub, sub, drop = FALSE])
      d <- dependence_metrics(R_rand)
      data.frame(MeanOffCor = d$MeanOffCor,
                 MaxOffCor = d$MaxOffCor,
                 MeanESS = d$MeanESS)
    }))
  } else {
    rand_dep <- NULL
  }

  # Build result object
  result <- structure(list(
    call = match.call(),
    type = type,
    selected = sel$selected,
    candidates = candidates,
    size = size,
    distance_metrics = data.frame(
      MinPD = obs_metrics$MinPD,
      MeanPD = obs_metrics$MeanPD,
      MeanNND = obs_metrics$MeanNND,
      MaxPD = obs_metrics$MaxPD,
      row.names = "observed"
    ),
    dependence_metrics = data.frame(
      MeanOffCor = dep$MeanOffCor,
      MaxOffCor = dep$MaxOffCor,
      MeanESS = dep$MeanESS,
      row.names = "observed"
    ),
    random_metrics = rand_metrics,
    random_dependence = rand_dep,
    empirical_p = emp_p,
    tree = tree,
    dist_mat = dist_mat,
    settings = list(
      type = type,
      n_random = n_random,
      cov_model = cov_model,
      lambda = lambda,
      half_life = half_life,
      seed = seed
    )
  ), class = "phylosubset_result")

  result
}

#' @export
#' @title Print method for phylosubset_result
#' @description Prints a concise summary of the analysis result.
#' @param x An object of class \code{phylosubset_result}.
#' @param ... Additional arguments (ignored).
print.phylosubset_result <- function(x, ...) {
  cat("PhyloSubset Result\n")
  cat("=================\n")
  cat("Type:", x$type, "\n")
  cat("Candidates:", length(x$candidates), "\n")
  cat("Subset size:", x$size, "\n")
  cat("Selected species:", length(x$selected), "\n")
  cat("\nDistance Metrics:\n")
  print(round(x$distance_metrics, 4))
  cat("\nDependence Metrics:\n")
  print(round(x$dependence_metrics, 4))
  if (!is.null(x$empirical_p)) {
    cat("\nEmpirical P-values:\n")
    print(x$empirical_p[, c("Metric", "Observed", "P_value", "SES")], digits = 3)
  }
  invisible(x)
}

#' @export
#' @title Summary method for phylosubset_result
#' @description Provides a detailed summary of the analysis result.
#' @param object An object of class \code{phylosubset_result}.
#' @param ... Additional arguments (ignored).
summary.phylosubset_result <- function(object, ...) {
  cat("PhyloSubset Analysis Summary\n")
  cat("============================\n")
  cat("Call:\n")
  print(object$call)
  cat("\nSettings:\n")
  cat("  Type:", object$settings$type, "\n")
  cat("  Candidates:", length(object$candidates), "\n")
  cat("  Subset size:", object$size, "\n")
  cat("  Random replicates:", object$settings$n_random, "\n")
  cat("  Covariance model:", object$settings$cov_model, "\n")

  cat("\nSelected species (", length(object$selected), "):\n", sep = "")
  cat("  ", paste(head(object$selected, 10), collapse = ", "),
      if (length(object$selected) > 10) ", ...", "\n", sep = "")

  cat("\nDistance Metrics:\n")
  dist_df <- object$distance_metrics
  if (!is.null(object$random_metrics)) {
    dist_summary <- data.frame(
      Observed = as.numeric(dist_df[1, ]),
      Baseline_Mean = sapply(names(dist_df), function(m) mean(object$random_metrics[[m]], na.rm = TRUE)),
      Baseline_SD = sapply(names(dist_df), function(m) stats::sd(object$random_metrics[[m]], na.rm = TRUE)),
      row.names = names(dist_df)
    )
    print(round(dist_summary, 4))
  } else {
    print(round(dist_df, 4))
  }

  cat("\nDependence Metrics:\n")
  dep_df <- object$dependence_metrics
  if (!is.null(object$random_dependence)) {
    dep_summary <- data.frame(
      Observed = as.numeric(dep_df[1, ]),
      Baseline_Mean = sapply(names(dep_df), function(m) mean(object$random_dependence[[m]], na.rm = TRUE)),
      Baseline_SD = sapply(names(dep_df), function(m) stats::sd(object$random_dependence[[m]], na.rm = TRUE)),
      row.names = names(dep_df)
    )
    print(round(dep_summary, 4))
  } else {
    print(round(dep_df, 4))
  }

  cat("\nEmpirical P-values:\n")
  print(object$empirical_p[, c("Metric", "Observed", "P_value", "SES")], digits = 3)

  cat("\nInterpretation:\n")
  if (object$type == "dispersed") {
    cat("  Higher MinPD/MeanPD/MeanNND indicates greater phylogenetic dispersion.\n")
    cat("  Lower MeanOffCor/MaxOffCor and higher MeanESS indicate weaker dependence.\n")
  } else {
    cat("  Lower MeanPD/MeanNND/MaxPD indicates tighter phylogenetic clustering.\n")
    cat("  Higher MeanOffCor/MaxOffCor and lower MeanESS indicate stronger dependence.\n")
  }

  invisible(object)
}

#' @export
#' @title Convert phylosubset_result to data frame
#' @description Extracts key results as a flat data frame.
#' @param x An object of class \code{phylosubset_result}.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Additional arguments (ignored).
as.data.frame.phylosubset_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    type = x$type,
    n_candidates = length(x$candidates),
    subset_size = x$size,
    MinPD = x$distance_metrics$MinPD,
    MeanPD = x$distance_metrics$MeanPD,
    MeanNND = x$distance_metrics$MeanNND,
    MaxPD = x$distance_metrics$MaxPD,
    MeanOffCor = x$dependence_metrics$MeanOffCor,
    MaxOffCor = x$dependence_metrics$MaxOffCor,
    MeanESS = x$dependence_metrics$MeanESS,
    cov_model = x$settings$cov_model,
    stringsAsFactors = FALSE
  )
}
