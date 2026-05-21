# 02_run_figure1_sanity_checks.R
# Figure 1: Idealized 32-tip tree sanity checks
# Balanced tree dispersed, balanced tree clustered, ladder tree clustered
# Fixed: n_tips = 32, s = 8
#
# Usage: Rscript scripts/02_run_figure1_sanity_checks.R
setwd("/home/huangr/projects/2026_04_distance")
source("R/01_load_modules.R")
load_project_modules()

cat("=== 02_run_figure1_sanity_checks ===\n")

N_TIPS <- 32
S <- 8
N_NULL <- 1000

out_dir <- file.path(RESULTS_DIR, "figure1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Generate trees ----
cat("Generating 32-tip ultrametric trees...\n")

balanced_tree <- make_ultrametric_balanced_tree_32()
check_sanity_tree(balanced_tree, "balanced")

ladder_tree <- make_ultrametric_ladder_tree(n_tips = 32)
check_sanity_tree(ladder_tree, "ladder")

# ---- Run analyses ----
cat("Running analyses (s =", S, ")...\n")

# Balanced dispersed
cat("  Balanced tree dispersed...\n")
dist_bal <- create_distance_object(balanced_tree)
disp_bal <- run_dispersed_algorithm(dist_bal, S)

# Balanced clustered
cat("  Balanced tree clustered...\n")
clust_bal <- select_best_clustered_neighborhood(dist_bal, S)

# Ladder clustered
cat("  Ladder tree clustered...\n")
dist_lad <- create_distance_object(ladder_tree)
clust_lad <- select_best_clustered_neighborhood(dist_lad, S)

# Random baseline (on balanced tree, same for all)
set.seed(GLOBAL_SEED)
random_idx <- sample_random_subsets(dist_bal, S, N_NULL, replace = FALSE)
random_dist <- calc_multiple_subsets_metrics_extended(dist_bal$dist_mat, random_idx)

# ---- Distance metrics ----
disp_bal_metrics <- calc_subset_metrics_extended(dist_bal$dist_mat, disp_bal$final_subset)
clust_bal_metrics <- calc_subset_metrics_extended(dist_bal$dist_mat, clust_bal$final_subset)
clust_lad_metrics <- calc_subset_metrics_extended(dist_lad$dist_mat, clust_lad$final_subset)

# ---- Summary table ----
summary_list <- list()

for (case_name in c("Balanced_Dispersed", "Balanced_Clustered", "Ladder_Clustered")) {
  obs <- switch(case_name,
    "Balanced_Dispersed" = disp_bal_metrics,
    "Balanced_Clustered" = clust_bal_metrics,
    "Ladder_Clustered" = clust_lad_metrics
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
      Case = case_name,
      Metric = mname,
      Observed = obs_val,
      Null_Mean = mean(null_vals, na.rm = TRUE),
      Null_SD = sd(null_vals, na.rm = TRUE),
      SES = ses_val,
      P_value = p_val,
      Tail = tail,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(out_dir, "figure1_sanity_summary_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_sanity_summary_s8.csv\n")

# ---- Selected species table ----
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

# ---- Distance metrics table ----
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

# ---- Tree plots ----
cat("Saving tree plots...\n")

save_tree_plot_pdf(balanced_tree, disp_bal$final_subset_names,
                   file.path(out_dir, "Fig1a_balanced_tree_dispersed_subset_s8.pdf"),
                   main = "Balanced Tree - Dispersed Subset (s=8)",
                   highlight_col = "red")

save_tree_plot_pdf(balanced_tree, clust_bal$final_subset_names,
                   file.path(out_dir, "Fig1b_balanced_tree_clustered_subset_s8.pdf"),
                   main = "Balanced Tree - Clustered Subset (s=8)",
                   highlight_col = "blue")

save_tree_plot_pdf(ladder_tree, clust_lad$final_subset_names,
                   file.path(out_dir, "Fig1c_ladder_tree_clustered_subset_s8.pdf"),
                   main = "Ladder Tree - Clustered Subset (s=8)",
                   highlight_col = "blue")

cat("Done. Figure 1 results saved to:", out_dir, "\n")
