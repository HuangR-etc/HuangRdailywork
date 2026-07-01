# 18_run_eb_teacher_grid_observed_trends.R
# Teacher-requested EB rho grid: observed dependence trends only
#
# Writes a file with the same format as
# results/covariance_sensitivity_eb_extreme/
#   cricetidae_covariance_sensitivity_EB_extreme_observed_trends.csv
# but for the rho grid:
#   0, +/-1, +/-2, ..., +/-10
#
# Usage:
#   Rscript scripts/18_run_eb_teacher_grid_observed_trends.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}

source("R/01_load_modules.R")
load_project_modules()

cat("=== 18_run_eb_teacher_grid_observed_trends ===\n")

raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
out_dir <- file.path(RESULTS_DIR, "covariance_sensitivity_eb_extreme")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rho_values <- sort(unique(c(0, seq(-10, -1, by = 1), seq(1, 10, by = 1))))
cov_grid <- COV_SENSITIVITY_GRID

all_observed <- list()

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

  for (r in rho_values) {
    cat("  EB teacher grid (rho =", r, ")...\n")

    V_eb <- make_eb_covariance(pool_tree, r)
    disp_dep <- calc_dependence_from_V(V_eb, disp_names)
    clust_dep <- calc_dependence_from_V(V_eb, clust_names)

    all_observed[[length(all_observed) + 1]] <- rbind(
      data.frame(
        Analysis = pool_label,
        N = N_i,
        s = s_i,
        Subset_Type = "dispersed",
        rho = r,
        MeanOffCor = disp_dep$off_mean,
        MaxOffCor = disp_dep$rmax,
        MIESS = disp_dep$neff_mean,
        stringsAsFactors = FALSE
      ),
      data.frame(
        Analysis = pool_label,
        N = N_i,
        s = s_i,
        Subset_Type = "clustered",
        rho = r,
        MeanOffCor = clust_dep$off_mean,
        MaxOffCor = clust_dep$rmax,
        MIESS = clust_dep$neff_mean,
        stringsAsFactors = FALSE
      )
    )
  }
}

if (length(all_observed) == 0) {
  stop("No teacher-grid EB observed trends were generated.")
}

observed_df <- do.call(rbind, all_observed)
observed_df <- observed_df[order(observed_df$N, observed_df$s, observed_df$Subset_Type, observed_df$rho), ]

out_file <- file.path(
  out_dir,
  "cricetidae_covariance_sensitivity_EB_teacher_grid_observed_trends.csv"
)

write.csv(
  observed_df,
  out_file,
  row.names = FALSE
)

cat("Done. Teacher-grid EB observed trends saved to:\n")
cat("  ", out_file, "\n", sep = "")
cat("Observed trend rows:", nrow(observed_df), "\n")
