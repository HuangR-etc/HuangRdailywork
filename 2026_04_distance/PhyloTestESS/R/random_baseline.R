#' Random baseline helpers
NULL

#' @export
random_subsets <- function(candidates, size, n = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  replicate(n, sample(candidates, size, replace = FALSE), simplify = FALSE)
}

#' @export
random_baseline <- function(dist_mat, candidates, size, n = 1000, seed = NULL) {
  subsets <- random_subsets(candidates, size, n = n, seed = seed)
  subset_metrics_table(dist_mat, subsets)
}

#' @export
empirical_p_value <- function(observed, random_values, direction = c("greater", "less")) {
  direction <- match.arg(direction)
  n_null <- length(random_values)

  extreme <- if (direction == "greater") {
    sum(random_values >= observed, na.rm = TRUE)
  } else {
    sum(random_values <= observed, na.rm = TRUE)
  }

  (extreme + 1) / (n_null + 1)
}

#' @export
compare_to_random <- function(observed_metrics, random_metrics,
                              type = c("dispersed", "clustered")) {
  type <- match.arg(type)
  metrics <- c("MinPD", "MeanPD", "MeanNND", "MaxPD")

  do.call(rbind, lapply(metrics, function(metric_name) {
    obs <- observed_metrics[[metric_name]]
    null_vals <- random_metrics[[metric_name]]
    direction <- if (type == "dispersed") "greater" else "less"

    data.frame(
      Metric = metric_name,
      Observed = obs,
      Baseline_Mean = mean(null_vals, na.rm = TRUE),
      Baseline_SD = stats::sd(null_vals, na.rm = TRUE),
      SES = if (stats::sd(null_vals, na.rm = TRUE) > 0) {
        (obs - mean(null_vals, na.rm = TRUE)) / stats::sd(null_vals, na.rm = TRUE)
      } else {
        NA_real_
      },
      P_value = empirical_p_value(obs, null_vals, direction = direction),
      Direction = direction,
      stringsAsFactors = FALSE
    )
  }))
}
