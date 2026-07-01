# 17_run_eb_extreme_dependence_check.R
# Extreme-case EB dependence validation
#
# Reuses the existing Cricetidae raw case RDS files and recomputes only the
# three dependence diagnostics under extreme EB rho settings. This is a
# lightweight trend-check script and does not modify the manuscript-ready
# covariance sensitivity outputs.
#
# Usage:
#   Rscript scripts/17_run_eb_extreme_dependence_check.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}

source("R/01_load_modules.R")
load_project_modules()

cat("=== 17_run_eb_extreme_dependence_check ===\n")

raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
out_dir <- file.path(RESULTS_DIR, "covariance_sensitivity_eb_extreme")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

extreme_rates <- EB_EXTREME_RATE_VALUES
cov_grid <- COV_SENSITIVITY_GRID

all_summaries <- list()
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
  random_names <- result$random_names

  for (r in extreme_rates) {
    cat("  EB extreme (rho =", r, ")...\n")

    V_eb <- make_eb_covariance(pool_tree, r)
    disp_dep <- calc_dependence_from_V(V_eb, disp_names)
    clust_dep <- calc_dependence_from_V(V_eb, clust_names)
    random_dep <- calc_multiple_dependence_from_V(V_eb, random_names)

    eb_summary <- summarize_dependence_observed_vs_random(
      pool_label = pool_label,
      N = N_i,
      subset_size = s_i,
      disp_dep = disp_dep,
      clust_dep = clust_dep,
      random_dep_metrics = random_dep,
      covariance_model = "EB_extreme",
      covariance_param = r,
      covariance_param_label = sprintf("EB_rho_%+.1f", r)
    )
    all_summaries[[length(all_summaries) + 1]] <- eb_summary

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

if (length(all_summaries) == 0) {
  stop("No EB extreme summaries were generated.")
}

summary_df <- do.call(rbind, all_summaries)
observed_df <- do.call(rbind, all_observed)
observed_df <- observed_df[order(observed_df$N, observed_df$s, observed_df$Subset_Type, observed_df$rho), ]

write.csv(
  summary_df,
  file.path(out_dir, "cricetidae_covariance_sensitivity_EB_extreme_summary.csv"),
  row.names = FALSE
)

write.csv(
  observed_df,
  file.path(out_dir, "cricetidae_covariance_sensitivity_EB_extreme_observed_trends.csv"),
  row.names = FALSE
)

cat("Done. EB extreme dependence results saved to:", out_dir, "\n")
cat("Summary rows:", nrow(summary_df), "\n")
cat("Observed trend rows:", nrow(observed_df), "\n")
