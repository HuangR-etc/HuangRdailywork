# 12_prediction_metric_ess.R
# Prediction-metric-based empirical effective sample size
#
# This module implements a simulation-based approach to calibrate the
# uncertainty of predictive-performance metrics (RMSE, MAE, R^2) against
# an independent-sample benchmark. The result is a metric-specific
# empirical effective sample size (ESS) that complements the mean-based
# MeanESS diagnostic.
#
# Key design choices:
# - Only two lambda scenarios: 0 (independent) and 1 (full BM correlation)
# - Trait y and prediction error e are simulated independently but share
#   the same correlation structure under a given lambda
# - Predictive R^2 is allowed to be negative (not truncated to 0)
# - ESS is obtained by interpolation against the lambda=0 benchmark curve
# - If target uncertainty exceeds the n=4 benchmark, ESS is reported as <4
# - If target uncertainty is below the n=32 benchmark, ESS is reported as >32

# ============================================================
# 4.1 Construct lambda-transformed correlation matrix
# ============================================================

#' Apply lambda transformation to a BM correlation matrix
#'
#' Off-diagonal elements are multiplied by lambda; diagonal is forced to 1.
#' lambda = 0 gives identity matrix; lambda = 1 gives full BM correlation.
#'
#' @param R_bm BM correlation matrix (symmetric, diag = 1)
#' @param lambda Scalar between 0 and 1
#' @return Transformed correlation matrix
make_lambda_correlation <- function(R_bm, lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("lambda must be a single numeric value.")
  }
  if (lambda < 0 || lambda > 1) {
    stop("lambda must be between 0 and 1.")
  }
  
  R_lam <- R_bm
  offdiag <- row(R_lam) != col(R_lam)
  R_lam[offdiag] <- lambda * R_lam[offdiag]
  diag(R_lam) <- 1
  R_lam <- (R_lam + t(R_lam)) / 2
  
  R_lam
}

# ============================================================
# 4.2 Stable generation of correlated multivariate normal
# ============================================================

#' Simulate correlated multivariate normal samples
#'
#' Uses Cholesky decomposition with fallback to eigen decomposition
#' for numerical stability.
#'
#' @param R Correlation matrix (n x n)
#' @param n_sim Number of simulation replicates
#' @param sd Standard deviation scaling factor
#' @param seed Optional random seed
#' @return Matrix (n_sim x n) of simulated values
simulate_correlated_normal <- function(R, n_sim, sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  R <- as.matrix(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  
  n <- nrow(R)
  
  # Prefer Cholesky, fall back to eigen if needed
  L <- tryCatch({
    t(chol(R))
  }, error = function(e) {
    eig <- eigen(R, symmetric = TRUE)
    vals <- pmax(eig$values, 0)
    eig$vectors %*% diag(sqrt(vals), nrow = length(vals))
  })
  
  Z <- matrix(rnorm(n_sim * n), nrow = n_sim, ncol = n)
  X <- Z %*% t(L)
  X <- sd * X
  
  X
}

# ============================================================
# 4.3 Simulate true values and predictions
# ============================================================

#' Simulate a prediction task with phylogenetic structure
#'
#' Generates true trait values y and prediction errors e independently,
#' both sharing the same correlation structure R. Predictions are
#' y_hat = y + e.
#'
#' @param R Correlation matrix (n x n)
#' @param n_sim Number of simulation replicates
#' @param trait_sd Standard deviation of the true trait
#' @param error_sd Standard deviation of the prediction error
#' @param seed Optional random seed
#' @return List with y, y_hat, and error matrices (each n_sim x n)
simulate_prediction_task <- function(R,
                                     n_sim = PRED_ESS_N_SIM,
                                     trait_sd = PRED_ESS_TRAIT_SD,
                                     error_sd = PRED_ESS_ERROR_SD,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  y <- simulate_correlated_normal(
    R = R,
    n_sim = n_sim,
    sd = trait_sd,
    seed = NULL
  )
  
  e <- simulate_correlated_normal(
    R = R,
    n_sim = n_sim,
    sd = error_sd,
    seed = NULL
  )
  
  y_hat <- y + e
  
  list(
    y = y,
    y_hat = y_hat,
    error = e
  )
}

# ============================================================
# 4.4 Calculate prediction metrics
# ============================================================

#' Calculate RMSE, MAE, and predictive R^2 from simulated data
#'
#' R^2 is the predictive R^2 (1 - SSE/SST) and is allowed to be negative.
#' Negative values are meaningful and should not be truncated.
#'
#' @param y Matrix of true values (n_sim x n)
#' @param y_hat Matrix of predicted values (n_sim x n)
#' @return Data frame with columns RMSE, MAE, R2
calc_prediction_metrics_matrix <- function(y, y_hat) {
  if (!all(dim(y) == dim(y_hat))) {
    stop("y and y_hat must have the same dimensions.")
  }
  
  residual <- y_hat - y
  
  rmse <- sqrt(rowMeans(residual^2))
  mae <- rowMeans(abs(residual))
  
  sse <- rowSums((y - y_hat)^2)
  y_centered <- y - rowMeans(y)
  sst <- rowSums(y_centered^2)
  
  r2 <- 1 - sse / sst
  r2[!is.finite(r2)] <- NA_real_
  
  data.frame(
    RMSE = rmse,
    MAE = mae,
    R2 = r2
  )
}

# ============================================================
# 4.5 Summarize metric uncertainty
# ============================================================

#' Summarize the empirical distribution of prediction metrics
#'
#' Computes mean, SD, quantiles, and the 95% empirical interval width
#' (Q97.5 - Q2.5) for each metric.
#'
#' @param metric_df Data frame with columns RMSE, MAE, R2
#' @param subset_type Character label for the subset type
#' @param lambda Lambda value used in simulation
#' @param benchmark_n Nominal sample size (for benchmark curve)
#' @param condition Character label for the simulation condition
#' @return Data frame with one row per metric
summarize_metric_uncertainty <- function(metric_df,
                                         subset_type,
                                         lambda,
                                         benchmark_n = NA_integer_,
                                         condition = NA_character_) {
  metrics <- c("RMSE", "MAE", "R2")
  
  out <- lapply(metrics, function(m) {
    x <- metric_df[[m]]
    x <- x[is.finite(x)]
    
    data.frame(
      Subset_Type = subset_type,
      Lambda = lambda,
      Benchmark_N = benchmark_n,
      Condition = condition,
      Metric = m,
      Mean = mean(x),
      SD = sd(x),
      Q025 = as.numeric(quantile(x, 0.025)),
      Q500 = as.numeric(quantile(x, 0.500)),
      Q975 = as.numeric(quantile(x, 0.975)),
      Interval_Width_95 = as.numeric(quantile(x, 0.975) - quantile(x, 0.025)),
      N_Sim = length(x),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

# ============================================================
# 4.6 Run lambda = 1 target subset simulation
# ============================================================

#' Run prediction-metric simulation for a target subset under lambda = 1
#'
#' Extracts the subset from the full BM correlation matrix, applies the
#' lambda transformation, simulates the prediction task, and summarizes
#' the metric uncertainty.
#'
#' @param R_bm_full Full BM correlation matrix (named)
#' @param subset_names Character vector of tip names in the subset
#' @param subset_type Character label (e.g., "dispersed", "clustered")
#' @param lambda Lambda value (typically 1 for phylogenetic scenario)
#' @param n_sim Number of simulation replicates
#' @param trait_sd Standard deviation of the true trait
#' @param error_sd Standard deviation of the prediction error
#' @param seed Random seed
#' @return List with subset_type, subset_names, lambda, R, metrics, summary
run_prediction_metric_target_subset <- function(R_bm_full,
                                                subset_names,
                                                subset_type,
                                                lambda = 1,
                                                n_sim = PRED_ESS_N_SIM,
                                                trait_sd = PRED_ESS_TRAIT_SD,
                                                error_sd = PRED_ESS_ERROR_SD,
                                                seed = PRED_ESS_SEED) {
  R_sub_bm <- R_bm_full[subset_names, subset_names, drop = FALSE]
  R_sub_lam <- make_lambda_correlation(R_sub_bm, lambda)
  
  sim <- simulate_prediction_task(
    R = R_sub_lam,
    n_sim = n_sim,
    trait_sd = trait_sd,
    error_sd = error_sd,
    seed = seed
  )
  
  metrics <- calc_prediction_metrics_matrix(sim$y, sim$y_hat)
  
  summary <- summarize_metric_uncertainty(
    metric_df = metrics,
    subset_type = subset_type,
    lambda = lambda,
    benchmark_n = length(subset_names),
    condition = paste0("lambda_", lambda, "_target")
  )
  
  list(
    subset_type = subset_type,
    subset_names = subset_names,
    lambda = lambda,
    R = R_sub_lam,
    metrics = metrics,
    summary = summary
  )
}


#' Run prediction-metric simulation for a target subset from any correlation matrix
#'
#' @param R_full Full named correlation matrix
#' @param subset_names Character vector of tip names in the subset
#' @param subset_type Character label
#' @param condition Character condition label
#' @param covariance_model Covariance model label
#' @param covariance_param Numeric covariance parameter
#' @param covariance_param_label Human-readable covariance parameter label
#' @param N Candidate pool size metadata
#' @param s Subset size metadata
#' @param n_sim Number of simulation replicates
#' @param trait_sd Standard deviation of the true trait
#' @param error_sd Standard deviation of the prediction error
#' @param seed Random seed
#' @return List with subset metadata, R, metrics, and summary
run_prediction_metric_target_subset_from_R <- function(R_full,
                                                       subset_names,
                                                       subset_type,
                                                       condition,
                                                       covariance_model,
                                                       covariance_param,
                                                       covariance_param_label,
                                                       N = NA_integer_,
                                                       s = length(subset_names),
                                                       n_sim = PRED_ESS_N_SIM,
                                                       trait_sd = PRED_ESS_TRAIT_SD,
                                                       error_sd = PRED_ESS_ERROR_SD,
                                                       seed = PRED_ESS_SEED) {
  R_sub <- R_full[subset_names, subset_names, drop = FALSE]
  R_sub <- (R_sub + t(R_sub)) / 2
  diag(R_sub) <- 1
  
  sim <- simulate_prediction_task(
    R = R_sub,
    n_sim = n_sim,
    trait_sd = trait_sd,
    error_sd = error_sd,
    seed = seed
  )
  
  metrics <- calc_prediction_metrics_matrix(sim$y, sim$y_hat)
  
  summary <- summarize_metric_uncertainty(
    metric_df = metrics,
    subset_type = subset_type,
    lambda = NA_real_,
    benchmark_n = length(subset_names),
    condition = condition
  )
  
  summary$N <- N
  summary$s <- s
  summary$Covariance_Model <- covariance_model
  summary$Covariance_Param <- covariance_param
  summary$Covariance_Param_Label <- covariance_param_label
  
  list(
    subset_type = subset_type,
    subset_names = subset_names,
    condition = condition,
    covariance_model = covariance_model,
    covariance_param = covariance_param,
    covariance_param_label = covariance_param_label,
    R = R_sub,
    metrics = metrics,
    summary = summary
  )
}

# ============================================================
# 4.7 Build lambda = 0 independent benchmark curve
# ============================================================

#' Build the independent-sample benchmark curve
#'
#' For each n in n_values, simulates the prediction task under lambda = 0
#' (identity correlation matrix) and summarizes metric uncertainty.
#' The resulting curve maps sample size to 95% empirical interval width.
#'
#' @param n_values Integer vector of sample sizes (e.g., 4:32)
#' @param n_sim Number of simulation replicates per n
#' @param trait_sd Standard deviation of the true trait
#' @param error_sd Standard deviation of the prediction error
#' @param seed Base random seed (offset by n_i for each n)
#' @return Data frame with benchmark curve (one row per metric per n)
run_independent_benchmark_curve <- function(n_values = PRED_ESS_BENCHMARK_N,
                                            n_sim = PRED_ESS_N_SIM,
                                            trait_sd = PRED_ESS_TRAIT_SD,
                                            error_sd = PRED_ESS_ERROR_SD,
                                            seed = PRED_ESS_SEED) {
  out <- list()
  
  for (n_i in n_values) {
    R_i <- diag(n_i)
    
    sim <- simulate_prediction_task(
      R = R_i,
      n_sim = n_sim,
      trait_sd = trait_sd,
      error_sd = error_sd,
      seed = seed + n_i
    )
    
    metrics <- calc_prediction_metrics_matrix(sim$y, sim$y_hat)
    
    out[[as.character(n_i)]] <- summarize_metric_uncertainty(
      metric_df = metrics,
      subset_type = "independent_benchmark",
      lambda = 0,
      benchmark_n = n_i,
      condition = "lambda_0_independent_benchmark"
    )
  }
  
  do.call(rbind, out)
}

# ============================================================
# 4.8 Monotonize benchmark curve widths
# ============================================================

#' Apply isotonic regression to enforce monotonicity on benchmark widths
#'
#' The 95% empirical interval width should decrease monotonically as
#' sample size increases. This function uses isotonic regression to
#' enforce this property, which stabilizes ESS interpolation.
#'
#' @param benchmark_df Data frame from run_independent_benchmark_curve
#' @return Data frame with added column Interval_Width_95_Monotone
monotonize_benchmark_widths <- function(benchmark_df) {
  metrics <- unique(benchmark_df$Metric)
  out <- list()
  
  for (m in metrics) {
    df_m <- benchmark_df[benchmark_df$Metric == m, ]
    df_m <- df_m[order(df_m$Benchmark_N), ]
    
    # Width should decrease as n increases.
    # Use isotonic regression on -width to get a nondecreasing curve,
    # then transform back.
    iso <- stats::isoreg(df_m$Benchmark_N, -df_m$Interval_Width_95)
    df_m$Interval_Width_95_Monotone <- -iso$yf
    
    out[[m]] <- df_m
  }
  
  do.call(rbind, out)
}

# ============================================================
# 4.9 Match prediction-metric-based ESS
# ============================================================

#' Estimate prediction-metric-based effective sample size
#'
#' For each target subset and metric, compares the 95% empirical interval
#' width under lambda = 1 to the independent-sample benchmark curve
#' (lambda = 0). The ESS is the benchmark sample size whose interval
#' width most closely matches the target.
#'
#' If the target width exceeds the n = 4 benchmark width, ESS is reported
#' as <4 (no extrapolation). If the target width is below the n = 32
#' benchmark width, ESS is reported as >32.
#'
#' @param target_summary Data frame from summarize_metric_uncertainty
#'   for target subsets under lambda = 1
#' @param benchmark_summary Data frame from run_independent_benchmark_curve
#'   (optionally monotonized)
#' @param use_monotone Logical; if TRUE, use monotonized widths
#' @return Data frame with one row per target subset per metric
estimate_prediction_metric_ess <- function(target_summary,
                                           benchmark_summary,
                                           use_monotone = TRUE) {
  if (use_monotone && !"Interval_Width_95_Monotone" %in% names(benchmark_summary)) {
    benchmark_summary <- monotonize_benchmark_widths(benchmark_summary)
  }
  
  width_col <- if (use_monotone) "Interval_Width_95_Monotone" else "Interval_Width_95"
  
  out <- list()
  
  for (i in seq_len(nrow(target_summary))) {
    row_i <- target_summary[i, ]
    metric_i <- row_i$Metric
    w_target <- row_i$Interval_Width_95
    
    bench_i <- benchmark_summary[benchmark_summary$Metric == metric_i, ]
    bench_i <- bench_i[order(bench_i$Benchmark_N), ]
    
    n_values <- bench_i$Benchmark_N
    widths <- bench_i[[width_col]]
    
    max_width <- max(widths, na.rm = TRUE)  # usually n = 4
    min_width <- min(widths, na.rm = TRUE)  # usually n = 32
    
    if (w_target > max_width) {
      ess <- NA_real_
      ess_label <- paste0("<", min(n_values))
      match_status <- "below_benchmark_range"
      ess_status <- paste0("<", min(n_values))
    } else if (w_target < min_width) {
      ess <- NA_real_
      ess_label <- paste0(">", max(n_values))
      match_status <- "above_benchmark_range"
      ess_status <- paste0(">", max(n_values))
    } else {
      # widths decrease as n increases. Reverse for approx x increasing.
      approx_df <- data.frame(width = widths, n = n_values)
      approx_df <- approx_df[order(approx_df$width), ]
      
      # Remove duplicated width values if needed
      approx_df <- approx_df[!duplicated(approx_df$width), ]
      
      ess <- as.numeric(stats::approx(
        x = approx_df$width,
        y = approx_df$n,
        xout = w_target,
        rule = 1
      )$y)
      
      ess_label <- sprintf("%.2f", ess)
      match_status <- "interpolated"
      ess_status <- "interpolated"
    }
    
    base_row <- data.frame(
      Subset_Type = row_i$Subset_Type,
      Metric = metric_i,
      Lambda_Target = row_i$Lambda,
      Nominal_N = row_i$Benchmark_N,
      Target_Interval_Width_95 = w_target,
      Interval_Width_95 = w_target,
      Prediction_Metric_ESS = ess,
      Prediction_Metric_ESS_Label = ess_label,
      Benchmark_N_Min = min(n_values),
      Benchmark_N_Max = max(n_values),
      Match_Status = match_status,
      ESS_Status = ess_status,
      stringsAsFactors = FALSE
    )
    
    metadata_cols <- c("N", "s", "Covariance_Model", "Covariance_Param", "Covariance_Param_Label", "Condition")
    for (col_i in metadata_cols) {
      if (col_i %in% names(row_i)) {
        base_row[[col_i]] <- row_i[[col_i]]
      }
    }
    
    out[[i]] <- base_row
  }
  
  do.call(rbind, out)
}


# ============================================================
# 4.10 PIESS significance, display, and direct-shift helpers
# ============================================================

p_to_stars_prediction <- function(p) {
  if (length(p) == 0 || is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

format_prediction_ess_display <- function(value,
                                          label = NA_character_,
                                          status = NA_character_) {
  if (!is.na(label) && nzchar(label)) return(label)
  if (!is.na(status) && status %in% c("<4", ">32")) return(status)
  if (is.finite(value)) return(sprintf("%.2f", value))
  "NA"
}

score_prediction_metric_ess <- function(ess_value,
                                        ess_label = NA_character_,
                                        benchmark_n_min = min(PRED_ESS_BENCHMARK_N),
                                        benchmark_n_max = max(PRED_ESS_BENCHMARK_N)) {
  if (is.finite(ess_value)) return(ess_value)

  if (!is.na(ess_label)) {
    if (grepl("^<", ess_label)) return(benchmark_n_min - 1e-6)
    if (grepl("^>", ess_label)) return(benchmark_n_max + 1e-6)
  }

  NA_real_
}

summarize_piess_observed_vs_random <- function(ess_summary,
                                               observed_subset_types = c("dispersed", "clustered"),
                                               random_subset_type = "random") {
  required_cols <- c(
    "N", "s", "Subset_Type", "Metric",
    "Interval_Width_95", "Prediction_Metric_ESS",
    "Prediction_Metric_ESS_Label"
  )
  missing_cols <- setdiff(required_cols, names(ess_summary))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ess_summary: ", paste(missing_cols, collapse = ", "))
  }

  group_cols <- intersect(
    c("N", "s", "Covariance_Model", "Covariance_Param", "Covariance_Param_Label", "Metric"),
    names(ess_summary)
  )

  ess_summary$PIESS_Score <- mapply(
    score_prediction_metric_ess,
    ess_summary$Prediction_Metric_ESS,
    ess_summary$Prediction_Metric_ESS_Label
  )

  observed <- ess_summary[ess_summary$Subset_Type %in% observed_subset_types, , drop = FALSE]
  random <- ess_summary[ess_summary$Subset_Type == random_subset_type, , drop = FALSE]

  out <- list()

  for (i in seq_len(nrow(observed))) {
    obs <- observed[i, , drop = FALSE]

    hit <- random
    for (g in group_cols) {
      hit <- hit[hit[[g]] == obs[[g]][1] | (is.na(hit[[g]]) & is.na(obs[[g]][1])), , drop = FALSE]
    }

    rand_scores <- hit$PIESS_Score
    rand_scores <- rand_scores[is.finite(rand_scores)]
    obs_score <- obs$PIESS_Score[1]

    if (length(rand_scores) == 0 || !is.finite(obs_score)) {
      p_val <- NA_real_
    } else if (obs$Subset_Type[1] == "dispersed") {
      p_val <- (1 + sum(rand_scores >= obs_score)) / (1 + length(rand_scores))
    } else if (obs$Subset_Type[1] == "clustered") {
      p_val <- (1 + sum(rand_scores <= obs_score)) / (1 + length(rand_scores))
    } else {
      p_val <- NA_real_
    }

    row_out <- obs
    row_out$Random_N <- length(rand_scores)
    row_out$Random_PIESS_Mean <- ifelse(length(rand_scores) > 0, mean(rand_scores), NA_real_)
    row_out$Random_PIESS_SD <- ifelse(length(rand_scores) > 1, sd(rand_scores), NA_real_)
    row_out$Random_PIESS_Q025 <- ifelse(length(rand_scores) > 0, as.numeric(quantile(rand_scores, 0.025)), NA_real_)
    row_out$Random_PIESS_Q975 <- ifelse(length(rand_scores) > 0, as.numeric(quantile(rand_scores, 0.975)), NA_real_)
    row_out$P_value <- p_val
    row_out$Stars <- p_to_stars_prediction(p_val)
    status_i <- if ("ESS_Status" %in% names(row_out)) row_out$ESS_Status else NA_character_
    row_out$PIESS_Display <- mapply(
      format_prediction_ess_display,
      row_out$Prediction_Metric_ESS,
      row_out$Prediction_Metric_ESS_Label,
      status_i
    )
    row_out$PIESS_Display_Stars <- paste0(row_out$PIESS_Display, row_out$Stars)

    out[[length(out) + 1]] <- row_out
  }

  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

calc_prediction_metric_shift <- function(target_summary,
                                         reference_condition = "lambda0_independent_target",
                                         scenario_condition = NULL,
                                         subset_type_filter = NULL) {
  required_cols <- c("N", "s", "Subset_Type", "Metric", "Mean", "Condition")
  missing_cols <- setdiff(required_cols, names(target_summary))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in target_summary: ", paste(missing_cols, collapse = ", "))
  }

  df <- target_summary
  if (!is.null(subset_type_filter)) {
    df <- df[df$Subset_Type %in% subset_type_filter, , drop = FALSE]
  }

  ref <- df[df$Condition == reference_condition, , drop = FALSE]
  scn <- df[df$Condition != reference_condition, , drop = FALSE]
  if (!is.null(scenario_condition)) {
    scn <- scn[scn$Condition %in% scenario_condition, , drop = FALSE]
  }

  group_cols <- intersect(c("N", "s", "Subset_Type", "Metric"), names(df))
  meta_cols <- intersect(
    c("Covariance_Model", "Covariance_Param", "Covariance_Param_Label", "Condition"),
    names(scn)
  )

  out <- list()

  for (i in seq_len(nrow(scn))) {
    row_i <- scn[i, , drop = FALSE]
    ref_i <- ref

    for (g in group_cols) {
      ref_i <- ref_i[ref_i[[g]] == row_i[[g]][1], , drop = FALSE]
    }

    if (nrow(ref_i) == 0) next

    independent_mean <- ref_i$Mean[1]
    scenario_mean <- row_i$Mean[1]

    denom <- abs(independent_mean)
    if (!is.finite(denom) || denom == 0) {
      percent_change <- NA_real_
    } else {
      percent_change <- 100 * (scenario_mean - independent_mean) / denom
    }

    out[[length(out) + 1]] <- data.frame(
      N = row_i$N[1],
      s = row_i$s[1],
      Subset_Type = row_i$Subset_Type[1],
      Metric = row_i$Metric[1],
      Independent_Mean = independent_mean,
      Scenario_Mean = scenario_mean,
      Absolute_Change = scenario_mean - independent_mean,
      Percent_Change = percent_change,
      row_i[, meta_cols, drop = FALSE],
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}
