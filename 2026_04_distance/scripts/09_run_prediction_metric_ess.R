# 09_run_prediction_metric_ess.R
# Run prediction-metric-based empirical effective sample size analysis
#
# This script reads the Cricetidae C512/s64 raw case (which contains the
# selected dispersed and clustered subsets), then:
# 1. Computes BM correlation matrix from the pool tree
# 2. Runs lambda=1 prediction-metric simulation for dispersed and clustered subsets
# 3. Builds lambda=0 independent benchmark curve for n = 4:32
# 4. Matches target interval widths to benchmark curve to estimate ESS
# 5. Saves results as RDS and CSV files
#
# Usage: Rscript scripts/09_run_prediction_metric_ess.R
setwd("/home/huangr/projects/2026_04_distance")
source("R/01_load_modules.R")
load_project_modules()

cat("=== 09_run_prediction_metric_ess ===\n")

out_dir <- file.path(RESULTS_DIR, "prediction_metric_ess")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_case_file <- file.path(
  RESULTS_DIR,
  "sensitivity",
  "raw_cases",
  "Cricetidae_C512_s64.rds"
)

if (!file.exists(raw_case_file)) {
  stop("Missing raw case RDS: ", raw_case_file,
       "\nRun scripts/04_run_cricetidae_sensitivity_grid.R --N 512 --s 64 --overwrite first.")
}

cat("Reading raw case from:", raw_case_file, "\n")
case <- readRDS(raw_case_file)

pool_tree <- case$pool_tree
disp_names <- case$dispersed$final_subset_names
clust_names <- case$clustered$final_subset_names

cat("Dispersed subset size:", length(disp_names), "\n")
cat("Clustered subset size:", length(clust_names), "\n")

# ---- Compute BM correlation matrix ----
cat("Computing BM covariance and correlation matrices...\n")
V_bm <- make_bm_covariance(pool_tree)
R_bm <- cov2cor(V_bm)
R_bm <- (R_bm + t(R_bm)) / 2
diag(R_bm) <- 1

# ---- Run target subset simulations ----
cat("Running target subset simulations (lambda = 1)...\n")

set.seed(PRED_ESS_SEED)

target_results <- list()

target_results[["dispersed"]] <- run_prediction_metric_target_subset(
  R_bm_full = R_bm,
  subset_names = disp_names,
  subset_type = "dispersed",
  lambda = PRED_ESS_LAMBDA_PHYLO,
  n_sim = PRED_ESS_N_SIM,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = PRED_ESS_ERROR_SD,
  seed = PRED_ESS_SEED + 1
)

target_results[["clustered"]] <- run_prediction_metric_target_subset(
  R_bm_full = R_bm,
  subset_names = clust_names,
  subset_type = "clustered",
  lambda = PRED_ESS_LAMBDA_PHYLO,
  n_sim = PRED_ESS_N_SIM,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = PRED_ESS_ERROR_SD,
  seed = PRED_ESS_SEED + 2
)

target_summary <- do.call(rbind, lapply(target_results, `[[`, "summary"))

cat("Target summary:\n")
print(target_summary)

# ---- Build independent benchmark curve ----
cat("Building independent benchmark curve (lambda = 0, n = 4:32)...\n")

benchmark_summary <- run_independent_benchmark_curve(
  n_values = PRED_ESS_BENCHMARK_N,
  n_sim = PRED_ESS_N_SIM,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = PRED_ESS_ERROR_SD,
  seed = PRED_ESS_SEED + 1000
)

cat("Benchmark summary (first few rows):\n")
print(head(benchmark_summary))

# ---- Monotonize benchmark curve ----
cat("Applying isotonic monotonicity correction to benchmark curve...\n")

benchmark_summary_mono <- monotonize_benchmark_widths(benchmark_summary)

# ---- Estimate prediction-metric ESS ----
cat("Estimating prediction-metric-based ESS...\n")

ess_summary <- estimate_prediction_metric_ess(
  target_summary = target_summary,
  benchmark_summary = benchmark_summary_mono,
  use_monotone = TRUE
)

cat("ESS summary:\n")
print(ess_summary)

# ---- Save results ----
cat("Saving results...\n")

result <- list(
  case_id = "Cricetidae_C512_s64_prediction_metric_ess",
  pool_label = "Cricetidae_C512",
  N = 512,
  s = 64,
  n_sim = PRED_ESS_N_SIM,
  benchmark_n = PRED_ESS_BENCHMARK_N,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = PRED_ESS_ERROR_SD,
  target_results = target_results,
  target_summary = target_summary,
  benchmark_summary = benchmark_summary,
  benchmark_summary_mono = benchmark_summary_mono,
  ess_summary = ess_summary
)

saveRDS(result, file.path(out_dir, "prediction_metric_ess_C512_s64.rds"))

write.csv(target_summary,
          file.path(out_dir, "prediction_metric_target_uncertainty_C512_s64.csv"),
          row.names = FALSE)

write.csv(benchmark_summary_mono,
          file.path(out_dir, "prediction_metric_independent_benchmark_C512_s64.csv"),
          row.names = FALSE)

write.csv(ess_summary,
          file.path(out_dir, "prediction_metric_ess_summary_C512_s64.csv"),
          row.names = FALSE)

cat("Done. Results saved to:", out_dir, "\n")
