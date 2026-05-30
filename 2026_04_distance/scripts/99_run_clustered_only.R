# 99_run_clustered_only.R
# Run only the clustered subset extraction using the new multi-start greedy-plus-exchange method
# This script runs all analyses that involve clustered subset selection and saves results.
#
# Usage: Rscript scripts/99_run_clustered_only.R
#        Rscript scripts/99_run_clustered_only.R --overwrite
setwd("/home/huangr/projects/2026_04_distance")
source("R/01_load_modules.R")
load_project_modules()

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
OVERWRITE <- "--overwrite" %in% args

cat("=== 99_run_clustered_only ===\n")
cat("OVERWRITE =", OVERWRITE, "\n")
cat("CLUSTERED_MAX_EXCHANGE_ITERATIONS =", CLUSTERED_MAX_EXCHANGE_ITERATIONS, "\n")
cat("CLUSTERED_EXCHANGE_TOL =", CLUSTERED_EXCHANGE_TOL, "\n\n")

# ============================================================
# Step 0: Ensure empirical trees are prepared
# ============================================================
cat("Step 0: Checking/Preparing empirical trees...\n")
if (!file.exists(file.path(PROCESSED_DIR, "carnivora_tree.rds")) ||
    !file.exists(file.path(PROCESSED_DIR, "cricetidae_full_tree.rds")) ||
    !file.exists(file.path(PROCESSED_DIR, "cricetidae_nested_pools.rds"))) {
  cat("  Empirical trees not found. Running 01_prepare_empirical_trees.R...\n")
  source("scripts/01_prepare_empirical_trees.R")
} else {
  cat("  Empirical trees already exist. Loading...\n")
  carnivora_tree <- readRDS(file.path(PROCESSED_DIR, "carnivora_tree.rds"))
  cricetidae_tree <- readRDS(file.path(PROCESSED_DIR, "cricetidae_full_tree.rds"))
  nested_pools <- readRDS(file.path(PROCESSED_DIR, "cricetidae_nested_pools.rds"))
}
cat("  Done.\n\n")

# ============================================================
# Step 1: Figure 1 sanity checks (clustered subsets)
# ============================================================
cat("Step 1: Figure 1 sanity checks...\n")

N_TIPS <- 32
S <- 8
N_NULL <- 1000

out_dir <- file.path(RESULTS_DIR, "figure1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

balanced_tree <- make_ultrametric_balanced_tree_32()
ladder_tree <- make_ultrametric_ladder_tree(n_tips = 32)

# Balanced dispersed (unchanged)
cat("  Balanced tree dispersed...\n")
dist_bal <- create_distance_object(balanced_tree)
disp_bal <- run_dispersed_algorithm(dist_bal, S)

# Balanced clustered (NEW method)
cat("  Balanced tree clustered (greedy+exchange)...\n")
clust_bal <- select_clustered_greedy_exchange(
  dist_obj = dist_bal,
  subset_size = S,
  max_iterations = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
  tol = CLUSTERED_EXCHANGE_TOL
)

# Ladder clustered (NEW method)
cat("  Ladder tree clustered (greedy+exchange)...\n")
dist_lad <- create_distance_object(ladder_tree)
clust_lad <- select_clustered_greedy_exchange(
  dist_obj = dist_lad,
  subset_size = S,
  max_iterations = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
  tol = CLUSTERED_EXCHANGE_TOL
)

# Random baselines (separate for balanced and ladder)
set.seed(GLOBAL_SEED)
random_idx_bal <- sample_random_subsets(dist_bal, S, N_NULL, replace = FALSE)
random_dist_bal <- calc_multiple_subsets_metrics_extended(dist_bal$dist_mat, random_idx_bal)

set.seed(GLOBAL_SEED)
random_idx_lad <- sample_random_subsets(dist_lad, S, N_NULL, replace = FALSE)
random_dist_lad <- calc_multiple_subsets_metrics_extended(dist_lad$dist_mat, random_idx_lad)

# Distance metrics
disp_bal_metrics <- calc_subset_metrics_extended(dist_bal$dist_mat, disp_bal$final_subset)
clust_bal_metrics <- calc_subset_metrics_extended(dist_bal$dist_mat, clust_bal$final_subset)
clust_lad_metrics <- calc_subset_metrics_extended(dist_lad$dist_mat, clust_lad$final_subset)

# Summary table
summary_list <- list()
for (case_name in c("Balanced_Dispersed", "Balanced_Clustered", "Ladder_Clustered")) {
  obs <- switch(case_name,
    "Balanced_Dispersed" = disp_bal_metrics,
    "Balanced_Clustered" = clust_bal_metrics,
    "Ladder_Clustered" = clust_lad_metrics
  )
  random_dist <- switch(case_name,
    "Balanced_Dispersed" = random_dist_bal,
    "Balanced_Clustered" = random_dist_bal,
    "Ladder_Clustered" = random_dist_lad
  )
  for (mname in c("MinPD", "MeanPD", "MeanNND", "MaxPD")) {
    obs_val <- obs[[mname]]
    null_vals <- random_dist[[mname]]
    ses_val <- calc_ses(obs_val, null_vals)
    if (case_name == "Balanced_Dispersed") {
      p_val <- calc_p_high(obs_val, null_vals)
      tail <- "upper"
    } else {
      p_val <- calc_p_low(obs_val, null_vals)
      tail <- "lower"
    }
    summary_list[[length(summary_list) + 1]] <- data.frame(
      Case = case_name, Metric = mname,
      Observed = obs_val, Null_Mean = mean(null_vals, na.rm = TRUE),
      Null_SD = sd(null_vals, na.rm = TRUE), SES = ses_val,
      P_value = p_val, Tail = tail, stringsAsFactors = FALSE
    )
  }
}
summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(out_dir, "figure1_sanity_summary_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_sanity_summary_s8.csv\n")

# Selected species table
max_len <- max(length(disp_bal$final_subset_names),
               length(clust_bal$final_subset_names),
               length(clust_lad$final_subset_names))
species_df <- data.frame(
  Species_Index = seq_len(max_len),
  Balanced_Dispersed = c(disp_bal$final_subset_names, rep(NA, max_len - length(disp_bal$final_subset_names))),
  Balanced_Clustered = c(clust_bal$final_subset_names, rep(NA, max_len - length(clust_bal$final_subset_names))),
  Ladder_Clustered = c(clust_lad$final_subset_names, rep(NA, max_len - length(clust_lad$final_subset_names))),
  stringsAsFactors = FALSE
)
write.csv(species_df, file.path(out_dir, "figure1_selected_tips_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_selected_tips_s8.csv\n")

# Distance metrics table
metrics_df <- data.frame(
  Case = c("Balanced_Dispersed", "Balanced_Clustered", "Ladder_Clustered"),
  MinPD = c(disp_bal_metrics$MinPD, clust_bal_metrics$MinPD, clust_lad_metrics$MinPD),
  MeanPD = c(disp_bal_metrics$MeanPD, clust_bal_metrics$MeanPD, clust_lad_metrics$MeanPD),
  MeanNND = c(disp_bal_metrics$MeanNND, clust_bal_metrics$MeanNND, clust_lad_metrics$MeanNND),
  MaxPD = c(disp_bal_metrics$MaxPD, clust_bal_metrics$MaxPD, clust_lad_metrics$MaxPD),
  stringsAsFactors = FALSE
)
write.csv(metrics_df, file.path(out_dir, "figure1_distance_metrics_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_distance_metrics_s8.csv\n")

# Save raw clustered results for later use
saveRDS(clust_bal, file.path(out_dir, "figure1_balanced_clustered_raw.rds"))
saveRDS(clust_lad, file.path(out_dir, "figure1_ladder_clustered_raw.rds"))
cat("  Saved: figure1 clustered raw RDS\n")

cat("  Done.\n\n")

# ============================================================
# Step 2: Carnivora main analysis
# ============================================================

# ============================================================
# Step 3: Cricetidae sensitivity grid
# ============================================================
cat("Step 3: Cricetidae sensitivity grid...\n")

nested_pools <- readRDS(file.path(PROCESSED_DIR, "cricetidae_nested_pools.rds"))

grid <- expand.grid(
  N = CRICETIDAE_POOL_SIZES,
  s = SENSITIVITY_SUBSET_SIZES
)
grid <- grid[grid$s < grid$N, ]
grid <- grid[order(grid$N, grid$s), ]

raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
sel_dir <- file.path(RESULTS_DIR, "sensitivity", "selected_subsets")
sum_dir <- file.path(RESULTS_DIR, "sensitivity", "summaries")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

all_dist_summaries <- list()
all_dep_summaries <- list()
all_selected_species <- list()

for (row_idx in seq_len(nrow(grid))) {
  N_i <- grid$N[row_idx]
  s_i <- grid$s[row_idx]
  
  pool_name <- paste0("C", N_i)
  pool_label <- paste0("Cricetidae_", pool_name)
  case_seed <- GLOBAL_SEED + N_i * 1000 + s_i
  
  cat("  Running: N =", N_i, ", s =", s_i, ", seed =", case_seed, "...\n")
  
  pool_tree <- nested_pools$pool_trees[[pool_name]]
  
  result <- run_empirical_case(
    pool_tree = pool_tree,
    pool_label = pool_label,
    subset_size = s_i,
    n_null_reps = N_NULL_REPS_DEFAULT,
    seed = case_seed,
    out_dir = raw_dir,
    save_raw = TRUE,
    overwrite = OVERWRITE
  )
  
  all_dist_summaries[[length(all_dist_summaries) + 1]] <- result$distance_summary
  all_dep_summaries[[length(all_dep_summaries) + 1]] <- result$dependence_summary_BM
  all_selected_species[[length(all_selected_species) + 1]] <- result$selected_species
}

dist_combined <- do.call(rbind, all_dist_summaries)
write.csv(dist_combined,
          file.path(sum_dir, "cricetidae_sensitivity_distance_summary.csv"), row.names = FALSE)

dep_combined <- do.call(rbind, all_dep_summaries)
write.csv(dep_combined,
          file.path(sum_dir, "cricetidae_sensitivity_dependence_BM_summary.csv"), row.names = FALSE)

species_combined <- do.call(rbind, all_selected_species)
write.csv(species_combined,
          file.path(sel_dir, "cricetidae_sensitivity_selected_species.csv"), row.names = FALSE)

cat("  Saved: Cricetidae sensitivity summaries\n")
cat("  Done.\n\n")

# ============================================================
# Step 4: Covariance sensitivity (reads raw RDS from Step 3)
# ============================================================
cat("Step 4: Covariance sensitivity...\n")

cov_out_dir <- file.path(RESULTS_DIR, "covariance_sensitivity")
dir.create(cov_out_dir, recursive = TRUE, showWarnings = FALSE)

# Load raw cases for the covariance sensitivity grid
cov_grid <- COV_SENSITIVITY_GRID

all_cov_summaries <- list()

for (row_idx in seq_len(nrow(cov_grid))) {
  N_i <- cov_grid$N[row_idx]
  s_i <- cov_grid$s[row_idx]
  
  pool_name <- paste0("Cricetidae_C", N_i)
  case_id <- paste0(pool_name, "_s", s_i)
  rds_file <- file.path(raw_dir, paste0(case_id, ".rds"))
  
  if (!file.exists(rds_file)) {
    cat("  WARNING: Raw case not found:", rds_file, "\n")
    cat("  Run Step 3 first or use --overwrite to regenerate.\n")
    next
  }
  
  cat("  Loading raw case:", case_id, "\n")
  raw_result <- readRDS(rds_file)
  
  pool_tree <- raw_result$pool_tree
  disp_names <- raw_result$dispersed$final_subset_names
  clust_names <- raw_result$clustered$final_subset_names
  random_names <- raw_result$random_names
  
  # BM (already computed, just extract)
  cat("    BM dependence...\n")
  bm_summary <- raw_result$dependence_summary_BM
  
  # lambda_BM
  cat("    lambda_BM dependence...\n")
  for (lambda_val in LAMBDA_VALUES) {
    V_lambda <- make_lambda_bm_covariance(pool_tree, lambda_val)
    
    disp_dep <- calc_dependence_from_V(V_lambda, disp_names)
    clust_dep <- calc_dependence_from_V(V_lambda, clust_names)
    random_dep <- calc_multiple_dependence_from_V(V_lambda, random_names)
    
    dep_summary <- summarize_dependence_observed_vs_random(
      pool_label = case_id, N = N_i, subset_size = s_i,
      disp_dep = disp_dep, clust_dep = clust_dep,
      random_dep_metrics = random_dep,
      covariance_model = "lambda_BM",
      covariance_param = lambda_val,
      covariance_param_label = paste0("lambda=", lambda_val)
    )
    
    all_cov_summaries[[length(all_cov_summaries) + 1]] <- dep_summary
  }
  
  # OU
  cat("    OU dependence...\n")
  for (hlf in OU_HALF_LIFE_FRACS) {
    V_ou <- make_ou_covariance(pool_tree, half_life_frac = hlf)
    
    disp_dep <- calc_dependence_from_V(V_ou, disp_names)
    clust_dep <- calc_dependence_from_V(V_ou, clust_names)
    random_dep <- calc_multiple_dependence_from_V(V_ou, random_names)
    
    dep_summary <- summarize_dependence_observed_vs_random(
      pool_label = case_id, N = N_i, subset_size = s_i,
      disp_dep = disp_dep, clust_dep = clust_dep,
      random_dep_metrics = random_dep,
      covariance_model = "OU",
      covariance_param = hlf,
      covariance_param_label = paste0("half_life_frac=", hlf)
    )
    
    all_cov_summaries[[length(all_cov_summaries) + 1]] <- dep_summary
  }
}

if (length(all_cov_summaries) > 0) {
  cov_combined <- do.call(rbind, all_cov_summaries)
  write.csv(cov_combined,
            file.path(cov_out_dir, "covariance_sensitivity_summary.csv"), row.names = FALSE)
  cat("  Saved: covariance_sensitivity_summary.csv\n")
} else {
  cat("  WARNING: No covariance sensitivity results generated.\n")
}

cat("  Done.\n\n")

# ============================================================
# Summary
# ============================================================
cat("========================================\n")
cat("All clustered subset extractions complete!\n")
cat("========================================\n")
cat("Results saved to:\n")
cat("  Figure 1:        ", file.path(RESULTS_DIR, "figure1"), "\n")
cat("  Carnivora:       ", file.path(RESULTS_DIR, "carnivora"), "\n")
cat("  Cricetidae sens: ", file.path(RESULTS_DIR, "sensitivity"), "\n")
cat("  Covariance sens: ", file.path(RESULTS_DIR, "covariance_sensitivity"), "\n")
cat("\n")
cat("Clustered method: multi-start greedy forward selection + one-for-one exchange refinement\n")
cat("Algorithm tag: clustered_multistart_greedy_exchange_meanpd_meannnd_maxpd\n")
cat("Exchange max iterations:", CLUSTERED_MAX_EXCHANGE_ITERATIONS, "\n")
cat("Exchange tolerance:", CLUSTERED_EXCHANGE_TOL, "\n")
cat("\n")
cat("To regenerate tables and figures, run:\n")
cat("  Rscript scripts/06_make_tables_and_figures.R\n")
cat("  Rscript scripts/07_make_metric_tables.R\n")
cat("  Rscript scripts/08_diagnose_s64_random_maxoffcor_pairs.R\n")
