# 15_run_prediction_metric_ess_full_grid.R
# BM full-grid PIESS.
#
# Usage examples:
# Rscript scripts/15_run_prediction_metric_ess_full_grid.R --N 512 --s 8 --include_random --n_random 5 --n_sim 100 --overwrite
# Rscript scripts/15_run_prediction_metric_ess_full_grid.R --include_random --n_random 1000 --n_sim 10000 --overwrite
# Rscript scripts/15_run_prediction_metric_ess_full_grid.R --N 512 --s 64 --no_random --overwrite
args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}
source("R/01_load_modules.R")
load_project_modules()

cat("=== 15_run_prediction_metric_ess_full_grid ===\n")

parse_cli_args <- function(args) {
  opts <- list(
    N = NULL,
    s = NULL,
    n_random = PRED_ESS_RANDOM_N_DEFAULT,
    n_sim = PRED_ESS_N_SIM,
    error_sd = PRED_ESS_ERROR_SD,
    error_var = PRED_ESS_ERROR_SD^2,
    include_random = FALSE,
    overwrite = FALSE,
    output_suffix = "_error_var_0p1"
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--N") {
      i <- i + 1; opts$N <- as.integer(args[i])
    } else if (a == "--s") {
      i <- i + 1; opts$s <- as.integer(args[i])
    } else if (a == "--n_random") {
      i <- i + 1; opts$n_random <- as.integer(args[i])
    } else if (a == "--n_sim") {
      i <- i + 1; opts$n_sim <- as.integer(args[i])
    } else if (a == "--error_sd") {
      i <- i + 1
      opts$error_sd <- as.numeric(args[i])
      opts$error_var <- opts$error_sd^2
    } else if (a == "--error_var") {
      i <- i + 1
      opts$error_var <- as.numeric(args[i])
      opts$error_sd <- sqrt(opts$error_var)
    } else if (a == "--output_suffix") {
      i <- i + 1; opts$output_suffix <- args[i]
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
                               covariance_param_label, N, s, n_sim,
                               error_sd, seed, random_id = NA_integer_) {
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
    error_sd = error_sd,
    seed = seed
  )
  summary_i <- target_result$summary
  summary_i$Random_ID <- random_id
  rm(target_result)
  gc(verbose = FALSE)
  summary_i
}

opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))
out_dir <- file.path(RESULTS_DIR, paste0("prediction_metric_ess_full_grid", opts$output_suffix))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rds_file <- file.path(out_dir, "prediction_metric_ess_full_grid.rds")
if (file.exists(rds_file) && !opts$overwrite) {
  stop("Output exists: ", rds_file, "\nUse --overwrite to replace it.")
}

raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
grid <- expand.grid(N = c(32, 64, 128, 256, 512), s = c(8, 16, 32, 64))
grid <- grid[grid$s < grid$N, , drop = FALSE]
if (!is.null(opts$N)) grid <- grid[grid$N == opts$N, , drop = FALSE]
if (!is.null(opts$s)) grid <- grid[grid$s == opts$s, , drop = FALSE]
if (nrow(grid) == 0) stop("No feasible N-s combinations selected.")

cat("Settings: n_sim=", opts$n_sim,
    ", error_sd=", signif(opts$error_sd, 6),
    ", error_var=", signif(opts$error_var, 6),
    ", include_random=", opts$include_random,
    ", n_random=", opts$n_random,
    ", output_suffix=", opts$output_suffix, "\n", sep = "")

cat("Building independent benchmark curve (lambda = 0, n = 4:32)...\n")
benchmark_summary <- run_independent_benchmark_curve(
  n_values = PRED_ESS_BENCHMARK_N,
  n_sim = opts$n_sim,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = opts$error_sd,
  seed = PRED_ESS_SEED + 2000
)
benchmark_summary_mono <- monotonize_benchmark_widths(benchmark_summary)

all_target_summary <- list()

for (row_idx in seq_len(nrow(grid))) {
  N_i <- grid$N[row_idx]
  s_i <- grid$s[row_idx]
  case_id <- paste0("Cricetidae_C", N_i, "_s", s_i)
  rds_in <- file.path(raw_dir, paste0(case_id, ".rds"))
  cat("Loading:", case_id, "...\n")
  if (!file.exists(rds_in)) {
    warning("Missing raw case, skipping: ", rds_in)
    next
  }

  case <- readRDS(rds_in)
  V_bm <- make_bm_covariance(case$pool_tree)
  R_bm <- cov2cor(V_bm)
  R_bm <- (R_bm + t(R_bm)) / 2
  diag(R_bm) <- 1

  subset_map <- list(
    dispersed = case$dispersed$final_subset_names,
    clustered = case$clustered$final_subset_names
  )

  for (subset_type in names(subset_map)) {
    subset_names <- subset_map[[subset_type]]
    R_ind <- diag(length(subset_names))
    rownames(R_ind) <- subset_names
    colnames(R_ind) <- subset_names
    seed_base <- PRED_ESS_SEED + N_i * 1000 + s_i * 10 + ifelse(subset_type == "clustered", 50000, 0)

    cat("  selected", subset_type, "BM lambda=1...\n")
    all_target_summary[[length(all_target_summary) + 1]] <- run_target_summary(
      R_full = R_bm,
      subset_names = subset_names,
      subset_type = subset_type,
      condition = "BM_lambda1_target",
      covariance_model = "BM",
      covariance_param = 1,
      covariance_param_label = "lambda=1",
      N = N_i,
      s = s_i,
      n_sim = opts$n_sim,
      error_sd = opts$error_sd,
      seed = seed_base + 1
    )

    cat("  selected", subset_type, "lambda=0 independent...\n")
    all_target_summary[[length(all_target_summary) + 1]] <- run_target_summary(
      R_full = R_ind,
      subset_names = subset_names,
      subset_type = subset_type,
      condition = "lambda0_independent_target",
      covariance_model = "independent",
      covariance_param = 0,
      covariance_param_label = "lambda=0",
      N = N_i,
      s = s_i,
      n_sim = opts$n_sim,
      error_sd = opts$error_sd,
      seed = seed_base + 2
    )
  }

  if (isTRUE(opts$include_random)) {
    random_list <- normalize_random_names(case$random_names)
    n_take <- min(opts$n_random, length(random_list))
    if (n_take == 0) warning("No random subsets available for ", case_id)
    for (random_id in seq_len(n_take)) {
      if (random_id %% 25 == 1 || random_id == n_take) {
        cat("  random BM lambda=1:", random_id, "/", n_take, "\n")
      }
      all_target_summary[[length(all_target_summary) + 1]] <- run_target_summary(
        R_full = R_bm,
        subset_names = random_list[[random_id]],
        subset_type = "random",
        condition = "BM_lambda1_target",
        covariance_model = "BM",
        covariance_param = 1,
        covariance_param_label = "lambda=1",
        N = N_i,
        s = s_i,
        n_sim = opts$n_sim,
        error_sd = opts$error_sd,
        seed = PRED_ESS_SEED + N_i * 1000 + s_i * 10 + 100000 + random_id,
        random_id = random_id
      )
    }
  }
}

if (length(all_target_summary) == 0) stop("No target summaries were generated.")
target_summary <- do.call(rbind, all_target_summary)

cat("Estimating prediction-metric-based ESS...\n")
piess_target_summary <- target_summary[target_summary$Condition == "BM_lambda1_target", , drop = FALSE]
ess_summary_all <- estimate_prediction_metric_ess(
  target_summary = piess_target_summary,
  benchmark_summary = benchmark_summary_mono,
  use_monotone = TRUE
)
if ("Random_ID" %in% names(piess_target_summary) && nrow(piess_target_summary) == nrow(ess_summary_all)) {
  ess_summary_all$Random_ID <- piess_target_summary$Random_ID
}

ess_summary_observed <- ess_summary_all[ess_summary_all$Subset_Type %in% c("dispersed", "clustered"), , drop = FALSE]
ess_summary_random <- ess_summary_all[ess_summary_all$Subset_Type == "random", , drop = FALSE]
ess_summary_vs_random <- summarize_piess_observed_vs_random(ess_summary_all)

ord <- c("N", "s", "Subset_Type", "Covariance_Model", "Covariance_Param", "Metric")
ess_summary_all <- ess_summary_all[do.call(order, ess_summary_all[intersect(ord, names(ess_summary_all))]), ]
ess_summary_observed <- ess_summary_observed[do.call(order, ess_summary_observed[intersect(ord, names(ess_summary_observed))]), ]
ess_summary_random <- ess_summary_random[do.call(order, ess_summary_random[intersect(ord, names(ess_summary_random))]), ]
ess_summary_vs_random <- ess_summary_vs_random[do.call(order, ess_summary_vs_random[intersect(ord, names(ess_summary_vs_random))]), ]

result <- list(
  case_grid = grid,
  n_sim = opts$n_sim,
  n_random = opts$n_random,
  include_random = opts$include_random,
  benchmark_n = PRED_ESS_BENCHMARK_N,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = opts$error_sd,
  error_var = opts$error_var,
  output_suffix = opts$output_suffix,
  target_summary = target_summary,
  benchmark_summary = benchmark_summary,
  benchmark_summary_mono = benchmark_summary_mono,
  ess_summary_all = ess_summary_all,
  ess_summary_observed = ess_summary_observed,
  ess_summary_random = ess_summary_random,
  ess_summary_vs_random = ess_summary_vs_random
)

cat("Saving results...\n")
saveRDS(result, rds_file)
write.csv(target_summary, file.path(out_dir, "prediction_metric_ess_full_grid_target_summary.csv"), row.names = FALSE)
write.csv(benchmark_summary_mono, file.path(out_dir, "prediction_metric_ess_full_grid_benchmark_summary.csv"), row.names = FALSE)
write.csv(ess_summary_all, file.path(out_dir, "prediction_metric_ess_full_grid_ess_summary_all.csv"), row.names = FALSE)
write.csv(ess_summary_observed, file.path(out_dir, "prediction_metric_ess_full_grid_ess_summary_observed.csv"), row.names = FALSE)
write.csv(ess_summary_random, file.path(out_dir, "prediction_metric_ess_full_grid_ess_summary_random.csv"), row.names = FALSE)
write.csv(ess_summary_vs_random, file.path(out_dir, "prediction_metric_ess_full_grid_ess_summary_vs_random.csv"), row.names = FALSE)
cat("Done. Full-grid PIESS results saved to:", out_dir, "\n")
