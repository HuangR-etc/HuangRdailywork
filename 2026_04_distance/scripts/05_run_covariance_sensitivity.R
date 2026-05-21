# 05_run_covariance_sensitivity.R
# Covariance sensitivity analysis
# Reads 4 representative Cricetidae raw case RDS files (already computed by
# 04_run_cricetidae_sensitivity_grid.R), and re-computes dependence diagnostics
# under lambda-transformed BM and OU covariance models.
# Does NOT re-select subsets or re-sample random baselines.
#
# Usage: Rscript scripts/05_run_covariance_sensitivity.R
setwd("/home/huangr/projects/2026_04_distance")
source("R/01_load_modules.R")
load_project_modules()

cat("=== 05_run_covariance_sensitivity ===\n")

raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
out_dir <- file.path(RESULTS_DIR, "covariance_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Representative cases ----
cov_grid <- COV_SENSITIVITY_GRID

all_summaries <- list()

for (row_idx in seq_len(nrow(cov_grid))) {
  N_i <- cov_grid$N[row_idx]
  s_i <- cov_grid$s[row_idx]
  
  pool_name <- paste0("C", N_i)
  pool_label <- paste0("Cricetidae_", pool_name)
  case_id <- paste0(pool_label, "_s", s_i)
  rds_file <- file.path(raw_dir, paste0(case_id, ".rds"))
  
  cat("Loading:", case_id, "...\n")
  
  if (!file.exists(rds_file)) {
    cat("  WARNING: File not found, skipping:", rds_file, "\n")
    next
  }
  
  result <- readRDS(rds_file)
  
  pool_tree <- result$pool_tree
  disp_names <- result$dispersed$final_subset_names
  clust_names <- result$clustered$final_subset_names
  random_names <- result$random_names
  
  # ---- BM covariance ----
  cat("  BM...\n")
  V_bm <- make_bm_covariance(pool_tree)
  
  disp_dep <- calc_dependence_from_V(V_bm, disp_names)
  clust_dep <- calc_dependence_from_V(V_bm, clust_names)
  random_dep <- calc_multiple_dependence_from_V(V_bm, random_names)
  
  bm_summary <- summarize_dependence_observed_vs_random(
    pool_label = pool_label,
    N = N_i,
    subset_size = s_i,
    disp_dep = disp_dep,
    clust_dep = clust_dep,
    random_dep_metrics = random_dep,
    covariance_model = "BM",
    covariance_param = NA_real_,
    covariance_param_label = "BM"
  )
  all_summaries[[length(all_summaries) + 1]] <- bm_summary
  
  # ---- Lambda-transformed BM ----
  for (lam in LAMBDA_VALUES) {
    cat("  lambda_BM (lambda =", lam, ")...\n")
    V_lam <- lambda_transform_cov(V_bm, lam)
    
    disp_dep <- calc_dependence_from_V(V_lam, disp_names)
    clust_dep <- calc_dependence_from_V(V_lam, clust_names)
    random_dep <- calc_multiple_dependence_from_V(V_lam, random_names)
    
    lam_summary <- summarize_dependence_observed_vs_random(
      pool_label = pool_label,
      N = N_i,
      subset_size = s_i,
      disp_dep = disp_dep,
      clust_dep = clust_dep,
      random_dep_metrics = random_dep,
      covariance_model = "lambda_BM",
      covariance_param = lam,
      covariance_param_label = paste0("lambda=", lam)
    )
    all_summaries[[length(all_summaries) + 1]] <- lam_summary
  }
  
  # ---- OU covariance ----
  for (hlf in OU_HALF_LIFE_FRACS) {
    cat("  OU (half_life_frac =", hlf, ")...\n")
    ou_result <- make_ou_covariance_by_half_life_fraction(pool_tree, hlf)
    
    disp_dep <- calc_dependence_from_V(ou_result$V, disp_names)
    clust_dep <- calc_dependence_from_V(ou_result$V, clust_names)
    random_dep <- calc_multiple_dependence_from_V(ou_result$V, random_names)
    
    ou_summary <- summarize_dependence_observed_vs_random(
      pool_label = pool_label,
      N = N_i,
      subset_size = s_i,
      disp_dep = disp_dep,
      clust_dep = clust_dep,
      random_dep_metrics = random_dep,
      covariance_model = "OU",
      covariance_param = hlf,
      covariance_param_label = paste0("half_life_frac=", hlf)
    )
    all_summaries[[length(all_summaries) + 1]] <- ou_summary
  }
}

# ---- Save combined summary ----
cat("Saving combined covariance sensitivity summary...\n")
combined <- do.call(rbind, all_summaries)
write.csv(combined,
          file.path(out_dir, "cricetidae_covariance_sensitivity_summary.csv"),
          row.names = FALSE)

cat("Done. Covariance sensitivity results saved to:", out_dir, "\n")
cat("Total rows:", nrow(combined), "\n")
