setwd("/home/galileo-group/galileouser05/huangr/2026_04_distance")

source("R/01_load_modules.R")
load_project_modules()

normalize_random_names <- function(x) {
  if (is.list(x)) return(x)
  if (is.matrix(x)) {
    return(lapply(seq_len(nrow(x)), function(i) x[i, ]))
  }
  stop("Unsupported random_names format.")
}

run_target_summary <- function(R_full, subset_names, subset_type, condition,
                               covariance_model, covariance_param,
                               covariance_param_label, N, s, n_sim, seed,
                               random_id = NA_integer_) {
  target_result <- run_prediction_metric_target_subset_from_R(
    R_full = R_full,
    subset_names = subset_names,
    subset_type = subset_type,
    condition = condition,
    covariance_model = covariance_model,
    covariance_param = covariance_param,
    covariance_param_label = covariance_param_label,
    N = N,
    s = s,
    n_sim = n_sim,
    trait_sd = PRED_ESS_TRAIT_SD,
    error_sd = PRED_ESS_ERROR_SD,
    seed = seed
  )
  summary_i <- target_result$summary
  summary_i$Random_ID <- random_id
  rm(target_result)
  gc(verbose = FALSE)
  summary_i
}

make_scenario_list <- function(pool_tree) {
  V_bm <- make_bm_covariance(pool_tree)
  scenarios <- list()

  for (lam in LAMBDA_VALUES) {
    scenarios[[length(scenarios) + 1]] <- list(
      model = "lambda_BM",
      param = lam,
      label = sprintf("lambda=%.2f", lam),
      condition = if (lam == 0) "lambda0_independent_target" else "covariance_sensitivity_target",
      R = cov2cor(lambda_transform_cov(V_bm, lam))
    )
  }

  for (hlf in OU_HALF_LIFE_FRACS) {
    ou_result <- make_ou_covariance_by_half_life_fraction(pool_tree, hlf)
    scenarios[[length(scenarios) + 1]] <- list(
      model = "OU",
      param = hlf,
      label = sprintf("h/H=%.2f", hlf),
      condition = "covariance_sensitivity_target",
      R = cov2cor(ou_result$V)
    )
  }

  for (rho in EB_RATE_VALUES) {
    V_eb <- make_eb_covariance(pool_tree, rho)
    scenarios[[length(scenarios) + 1]] <- list(
      model = "EB",
      param = rho,
      label = sprintf("rho=%.1f", rho),
      condition = "covariance_sensitivity_target",
      R = cov2cor(V_eb)
    )
  }

  lapply(scenarios, function(x) {
    x$R <- (x$R + t(x$R)) / 2
    diag(x$R) <- 1
    x
  })
}

recompute_covariance_eb <- function() {
  cat("=== EB-only recompute: covariance_sensitivity ===\n")

  raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
  out_dir <- file.path(RESULTS_DIR, "covariance_sensitivity")
  summary_csv <- file.path(out_dir, "cricetidae_covariance_sensitivity_summary.csv")
  figure_csv <- file.path(RESULTS_DIR, "figures", "Cricetidae_covariance_sensitivity_summary.csv")

  existing <- read.csv(summary_csv, stringsAsFactors = FALSE, check.names = FALSE)
  existing <- existing[existing$Covariance_Model != "EB", , drop = FALSE]

  eb_rows <- list()

  for (row_idx in seq_len(nrow(COV_SENSITIVITY_GRID))) {
    N_i <- COV_SENSITIVITY_GRID$N[row_idx]
    s_i <- COV_SENSITIVITY_GRID$s[row_idx]
    pool_name <- paste0("C", N_i)
    pool_label <- paste0("Cricetidae_", pool_name)
    case_id <- paste0(pool_label, "_s", s_i)
    rds_file <- file.path(raw_dir, paste0(case_id, ".rds"))

    cat("  Case:", case_id, "\n")
    case <- readRDS(rds_file)

    disp_names <- case$dispersed$final_subset_names
    clust_names <- case$clustered$final_subset_names
    random_names <- case$random_names

    for (rho in EB_RATE_VALUES) {
      cat("    EB rho =", rho, "\n")
      V_eb <- make_eb_covariance(case$pool_tree, rho)

      disp_dep <- calc_dependence_from_V(V_eb, disp_names)
      clust_dep <- calc_dependence_from_V(V_eb, clust_names)
      random_dep <- calc_multiple_dependence_from_V(V_eb, random_names)

      eb_rows[[length(eb_rows) + 1]] <- summarize_dependence_observed_vs_random(
        pool_label = pool_label,
        N = N_i,
        subset_size = s_i,
        disp_dep = disp_dep,
        clust_dep = clust_dep,
        random_dep_metrics = random_dep,
        covariance_model = "EB",
        covariance_param = rho,
        covariance_param_label = sprintf("EB_rho_%.1f", rho)
      )
    }
  }

  eb_df <- do.call(rbind, eb_rows)
  combined <- rbind(existing, eb_df)
  ord <- c("Analysis", "Pool", "N", "s", "Covariance_Model", "Covariance_Param", "Subset_Type", "Metric")
  combined <- combined[do.call(order, combined[intersect(ord, names(combined))]), , drop = FALSE]

  write.csv(combined, summary_csv, row.names = FALSE)
  write.csv(combined, figure_csv, row.names = FALSE)
  cat("  Updated:", summary_csv, "\n")
  cat("  Updated:", figure_csv, "\n")
}

recompute_piess_eb <- function() {
  cat("=== EB-only recompute: prediction_metric_ess_sensitivity ===\n")

  raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
  out_dir <- file.path(RESULTS_DIR, "prediction_metric_ess_sensitivity")
  piess_rds <- file.path(out_dir, "prediction_metric_ess_sensitivity.rds")

  existing <- readRDS(piess_rds)
  benchmark_summary_mono <- existing$benchmark_summary_mono

  all_target <- list()

  for (row_idx in seq_len(nrow(COV_SENSITIVITY_GRID))) {
    N_i <- COV_SENSITIVITY_GRID$N[row_idx]
    s_i <- COV_SENSITIVITY_GRID$s[row_idx]
    case_id <- paste0("Cricetidae_C", N_i, "_s", s_i)
    rds_in <- file.path(raw_dir, paste0(case_id, ".rds"))

    cat("  Case:", case_id, "\n")
    case <- readRDS(rds_in)

    subset_map <- list(
      dispersed = case$dispersed$final_subset_names,
      clustered = case$clustered$final_subset_names
    )
    random_list <- normalize_random_names(case$random_names)

    scenarios <- make_scenario_list(case$pool_tree)
    eb_idx <- which(vapply(scenarios, function(x) identical(x$model, "EB"), logical(1)))

    for (scenario_idx in eb_idx) {
      scenario <- scenarios[[scenario_idx]]
      cat("    EB", scenario$label, "(scenario_idx =", scenario_idx, ")\n")

      for (subset_type in names(subset_map)) {
        seed_i <- PRED_ESS_SEED + N_i * 1000 + s_i * 10 + scenario_idx +
          ifelse(subset_type == "clustered", 50000, 0)
        all_target[[length(all_target) + 1]] <- run_target_summary(
          R_full = scenario$R,
          subset_names = subset_map[[subset_type]],
          subset_type = subset_type,
          condition = scenario$condition,
          covariance_model = scenario$model,
          covariance_param = scenario$param,
          covariance_param_label = scenario$label,
          N = N_i,
          s = s_i,
          n_sim = existing$n_sim,
          seed = seed_i
        )
      }

      if (isTRUE(existing$include_random)) {
        n_take <- min(existing$n_random, length(random_list))
        for (random_id in seq_len(n_take)) {
          if (random_id %% 100 == 1 || random_id == n_take) {
            cat("      random", random_id, "/", n_take, "\n")
          }
          all_target[[length(all_target) + 1]] <- run_target_summary(
            R_full = scenario$R,
            subset_names = random_list[[random_id]],
            subset_type = "random",
            condition = scenario$condition,
            covariance_model = scenario$model,
            covariance_param = scenario$param,
            covariance_param_label = scenario$label,
            N = N_i,
            s = s_i,
            n_sim = existing$n_sim,
            seed = PRED_ESS_SEED + N_i * 1000 + s_i * 10 + scenario_idx * 100000 + random_id,
            random_id = random_id
          )
        }
      }
    }
  }

  eb_target <- do.call(rbind, all_target)
  eb_ess_all <- estimate_prediction_metric_ess(
    target_summary = eb_target,
    benchmark_summary = benchmark_summary_mono,
    use_monotone = TRUE
  )
  if ("Random_ID" %in% names(eb_target) && nrow(eb_target) == nrow(eb_ess_all)) {
    eb_ess_all$Random_ID <- eb_target$Random_ID
  }
  eb_ess_observed <- eb_ess_all[eb_ess_all$Subset_Type %in% c("dispersed", "clustered"), , drop = FALSE]
  eb_ess_random <- eb_ess_all[eb_ess_all$Subset_Type == "random", , drop = FALSE]
  eb_ess_vs_random <- summarize_piess_observed_vs_random(eb_ess_all)

  keep_non_eb <- function(df) df[df$Covariance_Model != "EB", , drop = FALSE]

  target_combined <- rbind(keep_non_eb(existing$target_summary), eb_target)
  ess_all_combined <- rbind(keep_non_eb(existing$ess_summary_all), eb_ess_all)
  ess_observed_combined <- rbind(keep_non_eb(existing$ess_summary_observed), eb_ess_observed)
  ess_random_combined <- rbind(keep_non_eb(existing$ess_summary_random), eb_ess_random)
  ess_vs_random_combined <- rbind(keep_non_eb(existing$ess_summary_vs_random), eb_ess_vs_random)

  ord <- c("N", "s", "Subset_Type", "Covariance_Model", "Covariance_Param", "Metric")
  target_combined <- target_combined[do.call(order, target_combined[intersect(ord, names(target_combined))]), , drop = FALSE]
  ess_all_combined <- ess_all_combined[do.call(order, ess_all_combined[intersect(ord, names(ess_all_combined))]), , drop = FALSE]
  ess_observed_combined <- ess_observed_combined[do.call(order, ess_observed_combined[intersect(ord, names(ess_observed_combined))]), , drop = FALSE]
  ess_random_combined <- ess_random_combined[do.call(order, ess_random_combined[intersect(ord, names(ess_random_combined))]), , drop = FALSE]
  ess_vs_random_combined <- ess_vs_random_combined[do.call(order, ess_vs_random_combined[intersect(ord, names(ess_vs_random_combined))]), , drop = FALSE]

  write.csv(target_combined, file.path(out_dir, "prediction_metric_ess_sensitivity_target_summary.csv"), row.names = FALSE)
  write.csv(benchmark_summary_mono, file.path(out_dir, "prediction_metric_ess_sensitivity_benchmark_summary.csv"), row.names = FALSE)
  write.csv(ess_all_combined, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_all.csv"), row.names = FALSE)
  write.csv(ess_observed_combined, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_observed.csv"), row.names = FALSE)
  write.csv(ess_random_combined, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_random.csv"), row.names = FALSE)
  write.csv(ess_vs_random_combined, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_vs_random.csv"), row.names = FALSE)
  write.csv(ess_observed_combined, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary.csv"), row.names = FALSE)

  existing$target_summary <- target_combined
  existing$ess_summary_all <- ess_all_combined
  existing$ess_summary_observed <- ess_observed_combined
  existing$ess_summary_random <- ess_random_combined
  existing$ess_summary_vs_random <- ess_vs_random_combined

  tmp_rds <- paste0(piess_rds, ".tmp")
  if (file.exists(tmp_rds)) file.remove(tmp_rds)
  saveRDS(existing, tmp_rds, compress = FALSE)
  if (file.exists(piess_rds)) file.remove(piess_rds)
  ok <- file.rename(tmp_rds, piess_rds)
  if (!isTRUE(ok)) stop("Failed to rename temp RDS into place.")

  cat("  Updated PIESS sensitivity files in:", out_dir, "\n")
}

recompute_covariance_eb()
recompute_piess_eb()
source("scripts/13_make_prediction_metric_ess_sensitivity_tables.R")
source("scripts/07_make_metric_tables.R")
cat("Done.\n")
