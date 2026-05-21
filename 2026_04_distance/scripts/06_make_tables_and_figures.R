# 06_make_tables_and_figures.R
# Read summaries and produce manuscript tables and figures.
# Does NOT run any analysis - only reads existing CSV/RDS files.
#
# Usage: Rscript scripts/06_make_tables_and_figures.R
setwd("/home/huangr/projects/2026_04_distance")
source("R/01_load_modules.R")
load_project_modules()

cat("=== 06_make_tables_and_figures ===\n")

fig_dir <- file.path(RESULTS_DIR, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Figure 1: Sanity check summary table
# ============================================================
cat("Reading Figure 1 results...\n")
fig1_summary <- read.csv(file.path(RESULTS_DIR, "figure1", "figure1_sanity_summary_s8.csv"))
fig1_metrics <- read.csv(file.path(RESULTS_DIR, "figure1", "figure1_distance_metrics_s8.csv"))
fig1_species <- read.csv(file.path(RESULTS_DIR, "figure1", "figure1_selected_tips_s8.csv"))

cat("  Figure 1 summary rows:", nrow(fig1_summary), "\n")

# ============================================================
# 2. Carnivora main analysis summary
# ============================================================
cat("Reading Carnivora results...\n")
carn_dist <- read.csv(file.path(RESULTS_DIR, "carnivora", "carnivora_distance_summary_s20.csv"))
carn_dep <- read.csv(file.path(RESULTS_DIR, "carnivora", "carnivora_dependence_BM_summary_s20.csv"))
carn_species <- read.csv(file.path(RESULTS_DIR, "carnivora", "carnivora_selected_species_s20.csv"))

cat("  Carnivora distance summary rows:", nrow(carn_dist), "\n")
cat("  Carnivora dependence summary rows:", nrow(carn_dep), "\n")

# ============================================================
# 3. Cricetidae sensitivity summaries
# ============================================================
cat("Reading Cricetidae sensitivity results...\n")
sens_dist <- read.csv(file.path(RESULTS_DIR, "sensitivity", "summaries",
                                 "cricetidae_sensitivity_distance_summary.csv"))
sens_dep <- read.csv(file.path(RESULTS_DIR, "sensitivity", "summaries",
                                "cricetidae_sensitivity_dependence_BM_summary.csv"))
sens_species <- read.csv(file.path(RESULTS_DIR, "sensitivity", "selected_subsets",
                                    "cricetidae_sensitivity_selected_species.csv"))

cat("  Sensitivity distance summary rows:", nrow(sens_dist), "\n")
cat("  Sensitivity dependence summary rows:", nrow(sens_dep), "\n")

# ============================================================
# 4. Covariance sensitivity summary
# ============================================================
cat("Reading covariance sensitivity results...\n")
cov_sens <- read.csv(file.path(RESULTS_DIR, "covariance_sensitivity",
                                "cricetidae_covariance_sensitivity_summary.csv"))
cat("  Covariance sensitivity rows:", nrow(cov_sens), "\n")

# ============================================================
# 5. Produce summary tables for manuscript
# ============================================================
cat("Producing manuscript summary tables...\n")

# Table 1: Carnivora distance + dependence combined
table1 <- merge(
  carn_dist[, c("Subset_Type", "Metric", "Observed", "SES", "P_value")],
  carn_dep[, c("Subset_Type", "Metric", "Observed", "SES", "P_value")],
  by = c("Subset_Type", "Metric"),
  all = TRUE,
  suffixes = c("_dist", "_dep")
)
write.csv(table1, file.path(fig_dir, "Table1_carnivora_summary.csv"), row.names = FALSE)
cat("  Saved: Table1_carnivora_summary.csv\n")

# Table 2: Cricetidae sensitivity distance (wide format)
# Pivot: N x s x Subset_Type -> columns for each Metric
table2 <- sens_dist[, c("N", "s", "Subset_Type", "Metric", "Observed", "SES", "P_value")]
write.csv(table2, file.path(fig_dir, "Table2_cricetidae_sensitivity_distance.csv"), row.names = FALSE)
cat("  Saved: Table2_cricetidae_sensitivity_distance.csv\n")

# Table 3: Cricetidae sensitivity dependence (wide format)
table3 <- sens_dep[, c("N", "s", "Subset_Type", "Metric", "Observed", "SES", "P_value")]
write.csv(table3, file.path(fig_dir, "Table3_cricetidae_sensitivity_dependence.csv"), row.names = FALSE)
cat("  Saved: Table3_cricetidae_sensitivity_dependence.csv\n")

# Table 4: Covariance sensitivity (neff_mean only, for compactness)
table4 <- cov_sens[cov_sens$Metric == "neff_mean",
                   c("N", "s", "Subset_Type", "Covariance_Model",
                     "Covariance_Param_Label", "Observed", "SES", "P_value")]
write.csv(table4, file.path(fig_dir, "Table4_covariance_sensitivity_neff.csv"), row.names = FALSE)
cat("  Saved: Table4_covariance_sensitivity_neff.csv\n")

# ============================================================
# 6. Summary statistics for manuscript text
# ============================================================
cat("Computing summary statistics for manuscript text...\n")

# Carnivora: dispersed vs random for MeanPD
carn_disp_meanpd <- carn_dist[carn_dist$Subset_Type == "dispersed" & carn_dist$Metric == "MeanPD", ]
cat(sprintf("Carnivora dispersed MeanPD: %.4f (SES = %.2f, p = %.4f)\n",
            carn_disp_meanpd$Observed, carn_disp_meanpd$SES, carn_disp_meanpd$P_value))

carn_clust_meanpd <- carn_dist[carn_dist$Subset_Type == "clustered" & carn_dist$Metric == "MeanPD", ]
cat(sprintf("Carnivora clustered MeanPD: %.4f (SES = %.2f, p = %.4f)\n",
            carn_clust_meanpd$Observed, carn_clust_meanpd$SES, carn_clust_meanpd$P_value))

# Carnivora dependence
carn_disp_neff <- carn_dep[carn_dep$Subset_Type == "dispersed" & carn_dep$Metric == "neff_mean", ]
cat(sprintf("Carnivora dispersed neff_mean: %.2f (SES = %.2f, p = %.4f)\n",
            carn_disp_neff$Observed, carn_disp_neff$SES, carn_disp_neff$P_value))

carn_clust_neff <- carn_dep[carn_dep$Subset_Type == "clustered" & carn_dep$Metric == "neff_mean", ]
cat(sprintf("Carnivora clustered neff_mean: %.2f (SES = %.2f, p = %.4f)\n",
            carn_clust_neff$Observed, carn_clust_neff$SES, carn_clust_neff$P_value))

cat("\nDone. All tables and figures prepared.\n")
cat("Output directory:", fig_dir, "\n")
