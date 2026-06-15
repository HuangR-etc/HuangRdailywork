# 02_run_figure1_sanity_checks.R
# Figure 1: Idealized 32-tip tree sanity checks
# Fixed: n_tips = 32, s = 8
#
# Usage: Rscript scripts/02_run_figure1_sanity_checks.R
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}
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

ladder_tree <- make_ultrametric_ladder_tree(n_tips = N_TIPS)
check_sanity_tree(ladder_tree, "ladder")

# ---- Run analyses ----
cat("Running analyses (s =", S, ")...\n")

dist_bal <- create_distance_object(balanced_tree)
dist_lad <- create_distance_object(ladder_tree)

cat("  Balanced tree dispersed...\n")
disp_bal <- run_dispersed_algorithm(dist_bal, S)

cat("  Balanced tree clustered...\n")
clust_bal <- select_clustered_greedy_exchange(
  dist_obj = dist_bal,
  subset_size = S,
  max_iterations = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
  tol = CLUSTERED_EXCHANGE_TOL
)

cat("  Ladder tree dispersed...\n")
disp_lad <- run_dispersed_algorithm(dist_lad, S)

cat("  Ladder tree clustered...\n")
clust_lad <- select_clustered_greedy_exchange(
  dist_obj = dist_lad,
  subset_size = S,
  max_iterations = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
  tol = CLUSTERED_EXCHANGE_TOL
)

# Random baseline (balanced and ladder trees separately)
set.seed(GLOBAL_SEED)
random_idx_bal <- sample_random_subsets(dist_bal, S, N_NULL, replace = FALSE)
random_dist_bal <- calc_multiple_subsets_metrics_extended(dist_bal$dist_mat, random_idx_bal)

set.seed(GLOBAL_SEED)
random_idx_lad <- sample_random_subsets(dist_lad, S, N_NULL, replace = FALSE)
random_dist_lad <- calc_multiple_subsets_metrics_extended(dist_lad$dist_mat, random_idx_lad)

# ---- Distance metrics ----
disp_bal_metrics <- calc_subset_metrics_extended(dist_bal$dist_mat, disp_bal$final_subset)
clust_bal_metrics <- calc_subset_metrics_extended(dist_bal$dist_mat, clust_bal$final_subset)
disp_lad_metrics <- calc_subset_metrics_extended(dist_lad$dist_mat, disp_lad$final_subset)
clust_lad_metrics <- calc_subset_metrics_extended(dist_lad$dist_mat, clust_lad$final_subset)

# ---- BM dependence diagnostics ----
V_bal <- make_bm_covariance(balanced_tree)
V_lad <- make_bm_covariance(ladder_tree)

dep_bal_disp <- calc_dependence_from_V(V_bal, disp_bal$final_subset_names)
dep_bal_clust <- calc_dependence_from_V(V_bal, clust_bal$final_subset_names)
dep_lad_disp <- calc_dependence_from_V(V_lad, disp_lad$final_subset_names)
dep_lad_clust <- calc_dependence_from_V(V_lad, clust_lad$final_subset_names)

cases <- list(
  Balanced_Dispersed = list(
    tree_type = "balanced", subset_type = "dispersed",
    metrics = disp_bal_metrics, dep = dep_bal_disp,
    names = disp_bal$final_subset_names, random = random_dist_bal
  ),
  Balanced_Clustered = list(
    tree_type = "balanced", subset_type = "clustered",
    metrics = clust_bal_metrics, dep = dep_bal_clust,
    names = clust_bal$final_subset_names, random = random_dist_bal
  ),
  Ladder_Dispersed = list(
    tree_type = "ladder", subset_type = "dispersed",
    metrics = disp_lad_metrics, dep = dep_lad_disp,
    names = disp_lad$final_subset_names, random = random_dist_lad
  ),
  Ladder_Clustered = list(
    tree_type = "ladder", subset_type = "clustered",
    metrics = clust_lad_metrics, dep = dep_lad_clust,
    names = clust_lad$final_subset_names, random = random_dist_lad
  )
)

# ---- Summary table ----
summary_list <- list()
for (case_name in names(cases)) {
  case_i <- cases[[case_name]]
  obs <- case_i$metrics
  random_dist <- case_i$random
  tail <- ifelse(case_i$subset_type == "dispersed", "upper", "lower")
  
  for (mname in c("MinPD", "MeanPD", "MeanNND", "MaxPD")) {
    obs_val <- obs[[mname]]
    null_vals <- random_dist[[mname]]
    p_val <- if (tail == "upper") calc_p_high(obs_val, null_vals) else calc_p_low(obs_val, null_vals)
    
    summary_list[[length(summary_list) + 1]] <- data.frame(
      Case = case_name,
      Tree_Type = case_i$tree_type,
      Subset_Type = case_i$subset_type,
      N = N_TIPS,
      s = S,
      Metric = mname,
      Observed = obs_val,
      Null_Mean = mean(null_vals, na.rm = TRUE),
      Null_SD = sd(null_vals, na.rm = TRUE),
      SES = calc_ses(obs_val, null_vals),
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
selected_long <- do.call(rbind, lapply(names(cases), function(case_name) {
  case_i <- cases[[case_name]]
  data.frame(
    Case = case_name,
    Tree_Type = case_i$tree_type,
    Subset_Type = case_i$subset_type,
    N = N_TIPS,
    s = S,
    Species_Index = seq_along(case_i$names),
    Species = case_i$names,
    stringsAsFactors = FALSE
  )
}))
write.csv(selected_long, file.path(out_dir, "figure1_selected_tips_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_selected_tips_s8.csv\n")

# ---- Distance metrics table ----
metrics_df <- do.call(rbind, lapply(names(cases), function(case_name) {
  case_i <- cases[[case_name]]
  data.frame(
    Case = case_name,
    Tree_Type = case_i$tree_type,
    Subset_Type = case_i$subset_type,
    N = N_TIPS,
    s = S,
    MinPD = case_i$metrics$MinPD,
    MeanPD = case_i$metrics$MeanPD,
    MeanNND = case_i$metrics$MeanNND,
    MaxPD = case_i$metrics$MaxPD,
    stringsAsFactors = FALSE
  )
}))
write.csv(metrics_df, file.path(out_dir, "figure1_distance_metrics_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_distance_metrics_s8.csv\n")

# ---- Dependence diagnostics table ----
dep_df <- do.call(rbind, lapply(names(cases), function(case_name) {
  case_i <- cases[[case_name]]
  data.frame(
    Case = case_name,
    Tree_Type = case_i$tree_type,
    Subset_Type = case_i$subset_type,
    N = N_TIPS,
    s = S,
    MeanOffCor = case_i$dep$off_mean,
    MaxOffCor = case_i$dep$rmax,
    MIESS = case_i$dep$neff_mean,
    off_mean = case_i$dep$off_mean,
    rmax = case_i$dep$rmax,
    MeanESS = case_i$dep$neff_mean,
    neff_mean = case_i$dep$neff_mean,
    stringsAsFactors = FALSE
  )
}))
write.csv(dep_df, file.path(out_dir, "figure1_dependence_BM_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_dependence_BM_s8.csv\n")

source_df <- merge(metrics_df, dep_df, by = c("Case", "Tree_Type", "Subset_Type", "N", "s"), all = TRUE)
write.csv(source_df, file.path(out_dir, "figure1_source_table_s8.csv"), row.names = FALSE)
cat("  Saved: figure1_source_table_s8.csv\n")

# ---- Tree plots ----
cat("Saving tree plots...\n")

panel_main <- function(case_name) {
  row_i <- dep_df[dep_df$Case == case_name, ]
  paste0(gsub("_", " - ", case_name), "\nMeanESS = ", sprintf("%.2f", row_i$MeanESS))
}

save_tree_plot_pdf(balanced_tree, disp_bal$final_subset_names,
                   file.path(out_dir, "Fig1a_balanced_tree_dispersed_subset_s8.pdf"),
                   main = panel_main("Balanced_Dispersed"),
                   highlight_col = "red")

save_tree_plot_pdf(balanced_tree, clust_bal$final_subset_names,
                   file.path(out_dir, "Fig1b_balanced_tree_clustered_subset_s8.pdf"),
                   main = panel_main("Balanced_Clustered"),
                   highlight_col = "blue")

save_tree_plot_pdf(ladder_tree, disp_lad$final_subset_names,
                   file.path(out_dir, "Fig1c_ladder_tree_dispersed_subset_s8.pdf"),
                   main = panel_main("Ladder_Dispersed"),
                   highlight_col = "red")

save_tree_plot_pdf(ladder_tree, clust_lad$final_subset_names,
                   file.path(out_dir, "Fig1d_ladder_tree_clustered_subset_s8.pdf"),
                   main = panel_main("Ladder_Clustered"),
                   highlight_col = "blue")

pdf(file.path(out_dir, "Figure1_idealized_tree_sanity_checks.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
plot_tree_with_subset(balanced_tree, disp_bal$final_subset_names,
                      main = panel_main("Balanced_Dispersed"), highlight_col = "red", cex = 0.45)
plot_tree_with_subset(balanced_tree, clust_bal$final_subset_names,
                      main = panel_main("Balanced_Clustered"), highlight_col = "blue", cex = 0.45)
plot_tree_with_subset(ladder_tree, disp_lad$final_subset_names,
                      main = panel_main("Ladder_Dispersed"), highlight_col = "red", cex = 0.45)
plot_tree_with_subset(ladder_tree, clust_lad$final_subset_names,
                      main = panel_main("Ladder_Clustered"), highlight_col = "blue", cex = 0.45)
dev.off()
cat("  Saved: Figure1_idealized_tree_sanity_checks.pdf\n")

cat("Done. Figure 1 results saved to:", out_dir, "\n")
