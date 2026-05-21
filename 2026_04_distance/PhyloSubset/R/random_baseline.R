#' Random baseline and empirical p-value functions
#'
#' Functions for generating random subsets, computing their metric
#' distributions, and calculating empirical p-values and standardized
#' effect sizes (SES).
#'
#' @name random_baseline
NULL

#' @rdname random_baseline
#' @export
#' @title Generate random subsets
#' @description Generates random subsets of a given size from a candidate pool.
#' @param candidates Character vector of candidate species names.
#' @param size Target subset size.
#' @param n Number of random subsets to generate.
#' @param seed Optional random seed for reproducibility.
#' @return A list of character vectors, each containing a random subset.
#' @examples
#' candidates <- paste0("sp", 1:20)
#' rand_subsets <- random_subsets(candidates, size = 5, n = 10)
#' length(rand_subsets)
random_subsets <- function(candidates, size, n = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  replicate(n, sample(candidates, size, replace = FALSE), simplify = FALSE)
}

#' @rdname random_baseline
#' @export
#' @title Compute random baseline metric distributions
#' @description Generates random subsets and computes distance metrics for
#'   each, returning the distribution of metrics.
#' @param dist_mat A symmetric patristic distance matrix.
#' @param candidates Character vector of candidate species names.
#' @param size Target subset size.
#' @param n Number of random subsets.
#' @param metrics Character vector of metrics to compute. Default is all four.
#' @param seed Optional random seed.
#' @return A data frame with columns for each metric and one row per random
#'   subset.
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(20)
#' D <- ape::cophenetic.phylo(tree)
#' baseline <- random_baseline(D, tree$tip.label, size = 5, n = 50)
#' head(baseline)
#' }
random_baseline <- function(dist_mat, candidates, size, n = 1000,
                            metrics = c("MinPD", "MeanPD", "MeanNND", "MaxPD"),
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  subsets <- random_subsets(candidates, size, n)

  results <- do.call(rbind, lapply(seq_len(n), function(i) {
    m <- distance_metrics(dist_mat, subsets[[i]])
    data.frame(
      MinPD = m$MinPD,
      MeanPD = m$MeanPD,
      MeanNND = m$MeanNND,
      MaxPD = m$MaxPD,
      stringsAsFactors = FALSE
    )
  }))

  results
}

#' @rdname random_baseline
#' @export
#' @title Calculate empirical p-value
#' @description Calculates a corrected empirical p-value as
#'   (b + 1) / (B + 1), where b is the number of random values as or more
#'   extreme than the observed value.
#' @param observed A single observed value.
#' @param random_values Numeric vector of random baseline values.
#' @param direction Direction of extremeness: \code{"greater"} (observed is
#'   larger than random) or \code{"less"} (observed is smaller than random).
#' @return A numeric p-value between 0 and 1.
#' @examples
#' set.seed(123)
#' obs <- 0.8
#' rand <- rnorm(100, mean = 0.5, sd = 0.2)
#' empirical_p_value(obs, rand, direction = "greater")
empirical_p_value <- function(observed, random_values,
                              direction = c("greater", "less")) {
  direction <- match.arg(direction)
  B <- length(random_values)

  if (direction == "greater") {
    b <- sum(random_values >= observed)
  } else {
    b <- sum(random_values <= observed)
  }

  (b + 1) / (B + 1)
}

#' @rdname random_baseline
#' @export
#' @title Compare observed subset to random baseline
#' @description Automatically compares observed distance metrics to a random
#'   baseline distribution, computing empirical p-values and standardized
#'   effect sizes (SES). The direction of comparison depends on the subset
#'   type: for dispersed, higher values are more extreme; for clustered,
#'   lower values are more extreme.
#' @param observed_metrics A list or data frame row with components MinPD,
#'   MeanPD, MeanNND, MaxPD.
#' @param random_metrics A data frame of random baseline metrics (as returned
#'   by \code{random_baseline}).
#' @param type Subset type: \code{"dispersed"} or \code{"clustered"}.
#' @return A data frame with columns Metric, Observed, Baseline_Mean,
#'   Baseline_SD, SES, P_value, and Direction.
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(20)
#' D <- ape::cophenetic.phylo(tree)
#' candidates <- tree$tip.label
#' obs <- distance_metrics(D, sample(candidates, 5))
#' rand <- random_baseline(D, candidates, size = 5, n = 50)
#' compare_to_random(obs, rand, type = "dispersed")
#' }
compare_to_random <- function(observed_metrics, random_metrics,
                              type = c("dispersed", "clustered")) {
  type <- match.arg(type)
  metrics_names <- c("MinPD", "MeanPD", "MeanNND", "MaxPD")

  results <- do.call(rbind, lapply(metrics_names, function(mname) {
    obs_val <- observed_metrics[[mname]]
    null_vals <- random_metrics[[mname]]

    ses <- if (stats::sd(null_vals) > 0) {
      (obs_val - mean(null_vals)) / stats::sd(null_vals)
    } else {
      NA_real_
    }

    if (type == "dispersed") {
      p_val <- empirical_p_value(obs_val, null_vals, "greater")
      dir <- "greater"
    } else {
      p_val <- empirical_p_value(obs_val, null_vals, "less")
      dir <- "less"
    }

    data.frame(
      Metric = mname,
      Observed = obs_val,
      Baseline_Mean = mean(null_vals, na.rm = TRUE),
      Baseline_SD = stats::sd(null_vals, na.rm = TRUE),
      SES = ses,
      P_value = p_val,
      Direction = dir,
      stringsAsFactors = FALSE
    )
  }))

  results
}
