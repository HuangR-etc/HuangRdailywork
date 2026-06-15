# 12_run_prediction_metric_ess_sensitivity.R
# Prediction-metric ESS sensitivity across representative N-s cases and covariance models.
#
# Usage examples:
# Rscript scripts/12_run_prediction_metric_ess_sensitivity.R --include_random --n_random 5 --n_sim 100 --overwrite
args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}
source("R/01_load_modules.R")
load_project_modules()

cat("=== 12_run_prediction_metric_ess_sensitivity ===\n")

parse_cli_args <- function(args) {
  opts <- list(
    n_random = PRED_ESS_RANDOM_N_DEFAULT,
    n_sim = PRED_ESS_N_SIM,
    include_random = FALSE,
    overwrite = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--n_random") {
      i <- i + 1; opts$n_random <- as.integer(args[i])
    } else if (a == "--n_sim") {
      i <- i + 1; opts$n_sim <- as.integer(args[i])
    } else if (a == "--include_random") {
      opts$include_random <- TRUE
    } else if (a == "--no_random") {
      opts$include_random <- FALSE
    } else if (a == "--overwrite") {
      opts$overwrite <- TRUE
    } else {
      stop("Unknown argument: ", a)
    }
    i <- i + 1
  }
  opts
}

normalize_random_names <- function(random_names) {
  if (is.null(random_names)) return(list())
  if (is.list(random_names) && !is.data.frame(random_names)) return(random_names)
  if (is.matrix(random_names) || is.data.frame(random_names)) {
    return(lapply(seq_len(nrow(random_names)), function(i) as.character(random_names[i, ])))
  }
  list(as.character(random_names))
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
  scenarios <- list(list(
    model = "BM",
    param = NA_real_,
    label = "BM",
    condition = "BM_target",
    R = cov2cor(V_bm)
  ))

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

  for (r in EB_RATE_VALUES) {
    V_eb <- make_eb_covariance_simple(pool_tree, r)
    scenarios[[length(scenarios) + 1]] <- list(
      model = "EB",
      param = r,
      label = sprintf("EB r=%.2f", r),
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

opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))
raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
out_dir <- file.path(RESULTS_DIR, "prediction_metric_ess_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rds_file <- file.path(out_dir, "prediction_metric_ess_sensitivity.rds")
if (file.exists(rds_file) && !opts$overwrite) {
  stop("Output exists: ", rds_file, "\nUse --overwrite to replace it.")
}

cov_grid <- COV_SENSITIVITY_GRID
cat("Settings: n_sim=", opts$n_sim,
    ", include_random=", opts$include_random,
    ", n_random=", opts$n_random, "\n", sep = "")

cat("Building independent benchmark curve (lambda = 0, n = 4:32)...\n")
benchmark_summary <- run_independent_benchmark_curve(
  n_values = PRED_ESS_BENCHMARK_N,
  n_sim = opts$n_sim,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = PRED_ESS_ERROR_SD,
  seed = PRED_ESS_SEED + 2000
)
benchmark_summary_mono <- monotonize_benchmark_widths(benchmark_summary)

all_target_summary <- list()

for (row_idx in seq_len(nrow(cov_grid))) {
  N_i <- cov_grid$N[row_idx]
  s_i <- cov_grid$s[row_idx]
  case_id <- paste0("Cricetidae_C", N_i, "_s", s_i)
  rds_in <- file.path(raw_dir, paste0(case_id, ".rds"))

  cat("Loading:", case_id, "...\n")
  if (!file.exists(rds_in)) {
    warning("File not found, skipping: ", rds_in)
    next
  }

  case <- readRDS(rds_in)
  subset_map <- list(
    dispersed = case$dispersed$final_subset_names,
    clustered = case$clustered$final_subset_names
  )
  scenarios <- make_scenario_list(case$pool_tree)

  for (scenario_idx in seq_along(scenarios)) {
    scenario <- scenarios[[scenario_idx]]
    cat("  ", scenario$model, scenario$label, "...\n")

    for (subset_type in names(subset_map)) {
      seed_i <- PRED_ESS_SEED + N_i * 1000 + s_i * 10 + scenario_idx + ifelse(subset_type == "clustered", 50000, 0)
      all_target_summary[[length(all_target_summary) + 1]] <- run_target_summary(
        R_full = scenario$R,
        subset_names = subset_map[[subset_type]],
        subset_type = subset_type,
        condition = scenario$condition,
        covariance_model = scenario$model,
        covariance_param = scenario$param,
        covariance_param_label = scenario$label,
        N = N_i,
        s = s_i,
        n_sim = opts$n_sim,
        seed = seed_i
      )
    }

    if (isTRUE(opts$include_random)) {
      random_list <- normalize_random_names(case$random_names)
      n_take <- min(opts$n_random, length(random_list))
      for (random_id in seq_len(n_take)) {
        if (random_id %% 25 == 1 || random_id == n_take) {
          cat("    random:", random_id, "/", n_take, "\n")
        }
        all_target_summary[[length(all_target_summary) + 1]] <- run_target_summary(
          R_full = scenario$R,
          subset_names = random_list[[random_id]],
          subset_type = "random",
          condition = scenario$condition,
          covariance_model = scenario$model,
          covariance_param = scenario$param,
          covariance_param_label = scenario$label,
          N = N_i,
          s = s_i,
          n_sim = opts$n_sim,
          seed = PRED_ESS_SEED + N_i * 1000 + s_i * 10 + scenario_idx * 100000 + random_id,
          random_id = random_id
        )
      }
    }
  }
}

if (length(all_target_summary) == 0) stop("No target summaries were generated.")
target_summary <- do.call(rbind, all_target_summary)

cat("Estimating prediction-metric-based ESS...\n")
ess_summary_all <- estimate_prediction_metric_ess(
  target_summary = target_summary,
  benchmark_summary = benchmark_summary_mono,
  use_monotone = TRUE
)
if ("Random_ID" %in% names(target_summary) && nrow(target_summary) == nrow(ess_summary_all)) {
  ess_summary_all$Random_ID <- target_summary$Random_ID
}

ess_summary_observed <- ess_summary_all[ess_summary_all$Subset_Type %in% c("dispersed", "clustered"), , drop = FALSE]
ess_summary_random <- ess_summary_all[ess_summary_all$Subset_Type == "random", , drop = FALSE]
ess_summary_vs_random <- summarize_piess_observed_vs_random(ess_summary_all)

cat("Calculating direct covariance metric shifts...\n")
metric_shift <- calc_prediction_metric_shift(
  target_summary = target_summary[target_summary$Subset_Type %in% c("dispersed", "clustered"), , drop = FALSE],
  reference_condition = "lambda0_independent_target"
)

# Include explicit zero rows for the lambda=0 reference, needed by panel tables.
ref_rows <- target_summary[
  target_summary$Subset_Type %in% c("dispersed", "clustered") &
    target_summary$Condition == "lambda0_independent_target",
  , drop = FALSE
]
if (nrow(ref_rows) > 0) {
  zero_shift <- data.frame(
    N = ref_rows$N,
    s = ref_rows$s,
    Subset_Type = ref_rows$Subset_Type,
    Metric = ref_rows$Metric,
    Independent_Mean = ref_rows$Mean,
    Scenario_Mean = ref_rows$Mean,
    Absolute_Change = 0,
    Percent_Change = 0,
    Covariance_Model = ref_rows$Covariance_Model,
    Covariance_Param = ref_rows$Covariance_Param,
    Covariance_Param_Label = ref_rows$Covariance_Param_Label,
    Condition = ref_rows$Condition,
    stringsAsFactors = FALSE
  )
  metric_shift <- rbind(zero_shift, metric_shift)
}

ord <- c("N", "s", "Subset_Type", "Covariance_Model", "Covariance_Param", "Metric")
ess_summary_all <- ess_summary_all[do.call(order, ess_summary_all[intersect(ord, names(ess_summary_all))]), ]
ess_summary_observed <- ess_summary_observed[do.call(order, ess_summary_observed[intersect(ord, names(ess_summary_observed))]), ]
ess_summary_random <- ess_summary_random[do.call(order, ess_summary_random[intersect(ord, names(ess_summary_random))]), ]
ess_summary_vs_random <- ess_summary_vs_random[do.call(order, ess_summary_vs_random[intersect(ord, names(ess_summary_vs_random))]), ]
target_summary <- target_summary[do.call(order, target_summary[intersect(ord, names(target_summary))]), ]

result <- list(
  case_grid = cov_grid,
  n_sim = opts$n_sim,
  n_random = opts$n_random,
  include_random = opts$include_random,
  benchmark_n = PRED_ESS_BENCHMARK_N,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = PRED_ESS_ERROR_SD,
  target_summary = target_summary,
  benchmark_summary = benchmark_summary,
  benchmark_summary_mono = benchmark_summary_mono,
  ess_summary_all = ess_summary_all,
  ess_summary_observed = ess_summary_observed,
  ess_summary_random = ess_summary_random,
  ess_summary_vs_random = ess_summary_vs_random,
  metric_shift = metric_shift
)

cat("Saving results...\n")
saveRDS(result, rds_file)
write.csv(target_summary, file.path(out_dir, "prediction_metric_ess_sensitivity_target_summary.csv"), row.names = FALSE)
write.csv(benchmark_summary_mono, file.path(out_dir, "prediction_metric_ess_sensitivity_benchmark_summary.csv"), row.names = FALSE)
write.csv(ess_summary_all, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_all.csv"), row.names = FALSE)
write.csv(ess_summary_observed, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_observed.csv"), row.names = FALSE)
write.csv(ess_summary_random, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_random.csv"), row.names = FALSE)
write.csv(ess_summary_vs_random, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_vs_random.csv"), row.names = FALSE)
write.csv(metric_shift, file.path(out_dir, "prediction_metric_direct_shift_covariance_sensitivity.csv"), row.names = FALSE)
write.csv(ess_summary_observed, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary.csv"), row.names = FALSE)
cat("Done. PIESS sensitivity results saved to:", out_dir, "\n")
