# 06_make_tables_and_figures.R
# Read existing Cricetidae summaries and produce unnumbered check tables.
# Does NOT run any analysis and does not depend on old Carnivora results.
#
# Usage: Rscript scripts/06_make_tables_and_figures.R
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}
source("R/01_load_modules.R")
load_project_modules()

cat("=== 06_make_tables_and_figures ===\n")

fig_dir <- file.path(RESULTS_DIR, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

# ============================================================
# 1. Figure 1 sanity-check sources
# ============================================================
cat("Reading Figure 1 results...\n")
fig1_summary <- read_required_csv(file.path(RESULTS_DIR, "figure1", "figure1_sanity_summary_s8.csv"))
fig1_metrics <- read_required_csv(file.path(RESULTS_DIR, "figure1", "figure1_distance_metrics_s8.csv"))
fig1_dep <- read_required_csv(file.path(RESULTS_DIR, "figure1", "figure1_dependence_BM_s8.csv"))
fig1_species <- read_required_csv(file.path(RESULTS_DIR, "figure1", "figure1_selected_tips_s8.csv"))

cat("  Figure 1 summary rows:", nrow(fig1_summary), "\n")
cat("  Figure 1 dependence rows:", nrow(fig1_dep), "\n")

# ============================================================
# 2. Cricetidae sensitivity summaries
# ============================================================
cat("Reading Cricetidae sensitivity results...\n")
sens_dist <- read_required_csv(file.path(
  RESULTS_DIR, "sensitivity", "summaries",
  "cricetidae_sensitivity_distance_summary.csv"
))
sens_dep <- read_required_csv(file.path(
  RESULTS_DIR, "sensitivity", "summaries",
  "cricetidae_sensitivity_dependence_BM_summary.csv"
))
sens_species <- read_required_csv(file.path(
  RESULTS_DIR, "sensitivity", "selected_subsets",
  "cricetidae_sensitivity_selected_species.csv"
))

cat("  Sensitivity distance summary rows:", nrow(sens_dist), "\n")
cat("  Sensitivity dependence summary rows:", nrow(sens_dep), "\n")
cat("  Selected species rows:", nrow(sens_species), "\n")

# ============================================================
# 3. Covariance sensitivity summary
# ============================================================
cat("Reading covariance sensitivity results...\n")
cov_sens <- read_required_csv(file.path(
  RESULTS_DIR, "covariance_sensitivity",
  "cricetidae_covariance_sensitivity_summary.csv"
))
cat("  Covariance sensitivity rows:", nrow(cov_sens), "\n")

# ============================================================
# 4. Produce unnumbered summary/check tables
# ============================================================
cat("Producing unnumbered summary/check tables...\n")

# Cricetidae BM dependence diagnostics.
bm_dep_summary <- sens_dep[, c("N", "s", "Subset_Type", "Metric", "Observed", "SES", "P_value")]
write.csv(bm_dep_summary, file.path(fig_dir, "Cricetidae_BM_dependence_summary.csv"), row.names = FALSE)
cat("  Saved: Cricetidae_BM_dependence_summary.csv\n")

# Covariance sensitivity, all requested dependence metrics.
cov_dep_summary <- cov_sens[cov_sens$Metric %in% c("off_mean", "rmax", "neff_mean"),
                            c("N", "s", "Subset_Type", "Covariance_Model",
                              "Covariance_Param_Label", "Metric", "Observed", "SES", "P_value")]
write.csv(cov_dep_summary, file.path(fig_dir, "Cricetidae_covariance_sensitivity_summary.csv"), row.names = FALSE)
cat("  Saved: Cricetidae_covariance_sensitivity_summary.csv\n")

# Figure 1 source table copy for manuscript assembly.
fig1_source <- merge(fig1_metrics, fig1_dep,
                     by = c("Case", "Tree_Type", "Subset_Type", "N", "s"),
                     all = TRUE)
write.csv(fig1_source, file.path(fig_dir, "Figure1_idealized_tree_sanity_checks_source.csv"), row.names = FALSE)
cat("  Saved: Figure1_idealized_tree_sanity_checks_source.csv\n")

cat("\nDone. Unnumbered summary/check tables prepared without Carnivora dependencies.\n")
cat("Output directory:", fig_dir, "\n")
