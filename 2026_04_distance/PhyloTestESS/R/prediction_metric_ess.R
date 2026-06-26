#' Prediction-metric-based ESS
#'
#' Simulation helpers for PIESS based on the 95% empirical interval width of
#' RMSE, MAE, and predictive R-squared.
NULL

PRED_ESS_N_SIM <- 10000
PRED_ESS_BENCHMARK_N <- 4:32
PRED_ESS_TRAIT_SD <- 1
PRED_ESS_ERROR_VAR <- 0.1
PRED_ESS_ERROR_SD <- sqrt(PRED_ESS_ERROR_VAR)
PRED_ESS_SEED <- 20260428 + 900000

make_lambda_correlation <- function(R_bm, lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0 || lambda > 1) {
    stop("lambda must be a single number between 0 and 1.", call. = FALSE)
  }

  R_lam <- R_bm
  offdiag <- row(R_lam) != col(R_lam)
  R_lam[offdiag] <- lambda * R_lam[offdiag]
  diag(R_lam) <- 1
  (R_lam + t(R_lam)) / 2
}

simulate_correlated_normal <- function(R, n_sim, sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  R <- as.matrix(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1

  L <- tryCatch(
    t(chol(R)),
    error = function(e) {
      eig <- eigen(R, symmetric = TRUE)
      eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))
    }
  )

  z <- matrix(rnorm(n_sim * nrow(R)), nrow = n_sim, ncol = nrow(R))
  sd * (z %*% t(L))
}

simulate_prediction_task <- function(R,
                                     n_sim = PRED_ESS_N_SIM,
                                     trait_sd = PRED_ESS_TRAIT_SD,
                                     error_sd = PRED_ESS_ERROR_SD,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  y <- simulate_correlated_normal(R = R, n_sim = n_sim, sd = trait_sd, seed = NULL)
  u <- simulate_correlated_normal(R = R, n_sim = n_sim, sd = error_sd, seed = NULL)

  list(
    y = y,
    y_hat = y + u,
    error = u
  )
}

calc_prediction_metrics_matrix <- function(y, y_hat) {
  residual <- y_hat - y
  rmse <- sqrt(rowMeans(residual^2))
  mae <- rowMeans(abs(residual))

  sse <- rowSums((y - y_hat)^2)
  sst <- rowSums((y - rowMeans(y))^2)
  r2 <- 1 - sse / sst
  r2[!is.finite(r2)] <- NA_real_

  data.frame(RMSE = rmse, MAE = mae, R2 = r2)
}

summarize_metric_uncertainty <- function(metric_df, subset_type,
                                         lambda,
                                         benchmark_n = NA_integer_,
                                         condition = NA_character_) {
  metrics <- c("RMSE", "MAE", "R2")

  do.call(rbind, lapply(metrics, function(metric_name) {
    x <- metric_df[[metric_name]]
    x <- x[is.finite(x)]

    data.frame(
      Subset_Type = subset_type,
      Lambda = lambda,
      Benchmark_N = benchmark_n,
      Condition = condition,
      Metric = metric_name,
      Mean = mean(x),
      SD = stats::sd(x),
      Q025 = as.numeric(stats::quantile(x, 0.025)),
      Q500 = as.numeric(stats::quantile(x, 0.5)),
      Q975 = as.numeric(stats::quantile(x, 0.975)),
      Interval_Width_95 = as.numeric(stats::quantile(x, 0.975) - stats::quantile(x, 0.025)),
      N_Sim = length(x),
      stringsAsFactors = FALSE
    )
  }))
}

#' @export
run_independent_benchmark_curve <- function(n_values = PRED_ESS_BENCHMARK_N,
                                            n_sim = PRED_ESS_N_SIM,
                                            trait_sd = PRED_ESS_TRAIT_SD,
                                            error_sd = PRED_ESS_ERROR_SD,
                                            seed = PRED_ESS_SEED) {
  out <- list()

  for (n_i in n_values) {
    sim <- simulate_prediction_task(
      R = diag(n_i),
      n_sim = n_sim,
      trait_sd = trait_sd,
      error_sd = error_sd,
      seed = seed + n_i
    )
    out[[as.character(n_i)]] <- summarize_metric_uncertainty(
      metric_df = calc_prediction_metrics_matrix(sim$y, sim$y_hat),
      subset_type = "independent_benchmark",
      lambda = 0,
      benchmark_n = n_i,
      condition = "lambda_0_independent_benchmark"
    )
  }

  do.call(rbind, out)
}

monotonize_benchmark_widths <- function(benchmark_df) {
  metrics <- unique(benchmark_df$Metric)

  do.call(rbind, lapply(metrics, function(metric_name) {
    df <- benchmark_df[benchmark_df$Metric == metric_name, , drop = FALSE]
    df <- df[order(df$Benchmark_N), , drop = FALSE]
    iso <- stats::isoreg(df$Benchmark_N, -df$Interval_Width_95)
    df$Interval_Width_95_Monotone <- -iso$yf
    df
  }))
}

#' @export
estimate_prediction_metric_ess <- function(target_summary, benchmark_summary,
                                           use_monotone = TRUE) {
  if (use_monotone && !"Interval_Width_95_Monotone" %in% names(benchmark_summary)) {
    benchmark_summary <- monotonize_benchmark_widths(benchmark_summary)
  }

  width_col <- if (use_monotone) "Interval_Width_95_Monotone" else "Interval_Width_95"

  do.call(rbind, lapply(seq_len(nrow(target_summary)), function(i) {
    target <- target_summary[i, , drop = FALSE]
    bench <- benchmark_summary[benchmark_summary$Metric == target$Metric[1], , drop = FALSE]
    bench <- bench[order(bench$Benchmark_N), , drop = FALSE]

    widths <- bench[[width_col]]
    n_values <- bench$Benchmark_N
    w_target <- target$Interval_Width_95[1]

    if (w_target > max(widths, na.rm = TRUE)) {
      ess <- NA_real_
      label <- paste0("<", min(n_values))
      status <- "below_benchmark_range"
      ess_status <- paste0("<", min(n_values))
    } else if (w_target < min(widths, na.rm = TRUE)) {
      ess <- NA_real_
      label <- paste0(">", max(n_values))
      status <- "above_benchmark_range"
      ess_status <- paste0(">", max(n_values))
    } else {
      approx_df <- data.frame(width = widths, n = n_values)
      approx_df <- approx_df[order(approx_df$width), , drop = FALSE]
      approx_df <- approx_df[!duplicated(approx_df$width), , drop = FALSE]
      ess <- as.numeric(stats::approx(
        x = approx_df$width,
        y = approx_df$n,
        xout = w_target,
        rule = 1
      )$y)
      label <- sprintf("%.2f", ess)
      status <- "interpolated"
      ess_status <- "interpolated"
    }

    out <- data.frame(
      Subset_Type = target$Subset_Type[1],
      Metric = target$Metric[1],
      Lambda_Target = target$Lambda[1],
      Nominal_N = target$Benchmark_N[1],
      Target_Interval_Width_95 = w_target,
      Interval_Width_95 = w_target,
      Prediction_Metric_ESS = ess,
      Prediction_Metric_ESS_Label = label,
      Benchmark_N_Min = min(n_values),
      Benchmark_N_Max = max(n_values),
      Match_Status = status,
      ESS_Status = ess_status,
      stringsAsFactors = FALSE
    )

    metadata_cols <- c("N", "s", "Covariance_Model", "Covariance_Param", "Covariance_Param_Label", "Condition")
    for (col_i in metadata_cols) {
      if (col_i %in% names(target)) {
        out[[col_i]] <- target[[col_i]]
      }
    }
    out
  }))
}

#' @export
analyze_prediction_metric_ess <- function(R_full, subset_names,
                                          subset_type = "target",
                                          n_sim = PRED_ESS_N_SIM,
                                          trait_sd = PRED_ESS_TRAIT_SD,
                                          error_sd = PRED_ESS_ERROR_SD,
                                          benchmark_n = NULL,
                                          seed = PRED_ESS_SEED,
                                          benchmark_seed = PRED_ESS_SEED + 2000,
                                          condition = "phylogenetic_target") {
  R_sub <- R_full[subset_names, subset_names, drop = FALSE]
  R_sub <- (R_sub + t(R_sub)) / 2
  diag(R_sub) <- 1

  if (is.null(benchmark_n)) {
    benchmark_n <- PRED_ESS_BENCHMARK_N
  }

  sim <- simulate_prediction_task(
    R = R_sub,
    n_sim = n_sim,
    trait_sd = trait_sd,
    error_sd = error_sd,
    seed = seed
  )

  target_summary <- summarize_metric_uncertainty(
    metric_df = calc_prediction_metrics_matrix(sim$y, sim$y_hat),
    subset_type = subset_type,
    lambda = NA_real_,
    benchmark_n = length(subset_names),
    condition = condition
  )

  benchmark_summary <- run_independent_benchmark_curve(
    n_values = benchmark_n,
    n_sim = n_sim,
    trait_sd = trait_sd,
    error_sd = error_sd,
    seed = benchmark_seed
  )

  list(
    subset_names = subset_names,
    target_summary = target_summary,
    benchmark_summary = benchmark_summary,
    piess = estimate_prediction_metric_ess(target_summary, benchmark_summary)
  )
}
