args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}

source("R/01_load_modules.R")
load_project_modules()

cat("=== 16_remap_prediction_metric_ess_sensitivity_benchmark64 ===\n")

out_dir <- file.path(RESULTS_DIR, "prediction_metric_ess_sensitivity")
target_file <- file.path(out_dir, "prediction_metric_ess_sensitivity_target_summary.csv")
rds_file <- file.path(out_dir, "prediction_metric_ess_sensitivity.rds")
vs_random_file <- file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_vs_random.csv")

if (!file.exists(target_file)) {
  stop("Missing target summary: ", target_file)
}
if (!file.exists(vs_random_file)) {
  stop("Missing observed-vs-random summary: ", vs_random_file)
}

cat("Reading target summary from:", target_file, "\n")
target_summary <- read.csv(target_file, stringsAsFactors = FALSE)
old_vs_random <- read.csv(vs_random_file, stringsAsFactors = FALSE)

affected_groups <- unique(
  old_vs_random[
    old_vs_random$Prediction_Metric_ESS_Label %in% c(">32", ">64"),
    c("N", "s", "Metric", "Covariance_Model", "Covariance_Param")
  ]
)

if (nrow(affected_groups) == 0) {
  cat("No observed rows are currently labeled >32 or >64. Nothing to remap.\n")
  quit(save = "no", status = 0)
}

cat(
  "Affected observed rows:",
  sum(old_vs_random$Prediction_Metric_ESS_Label %in% c(">32", ">64")),
  "\n"
)
cat("Affected groups:", nrow(affected_groups), "\n")

key_cols <- c("N", "s", "Metric", "Covariance_Model", "Covariance_Param")
make_key <- function(df) {
  do.call(
    paste,
    c(
      lapply(key_cols, function(col) {
        if (col %in% names(df)) {
          if (is.numeric(df[[col]])) sprintf("%.10f", df[[col]]) else as.character(df[[col]])
        } else {
          rep("", nrow(df))
        }
      }),
      sep = "|"
    )
  )
}

affected_keys <- make_key(affected_groups)
target_keys <- make_key(target_summary)

affected_target_summary <- target_summary[target_keys %in% affected_keys, , drop = FALSE]
unaffected_target_summary <- target_summary[!(target_keys %in% affected_keys), , drop = FALSE]

cat("Rows to remap in target_summary:", nrow(affected_target_summary), "\n")
cat("Rows left untouched in target_summary:", nrow(unaffected_target_summary), "\n")

benchmark_n_values <- 4:70
cat("Building extended independent benchmark curve (lambda = 0, n = 4:70)...\n")
benchmark_summary <- run_independent_benchmark_curve(
  n_values = benchmark_n_values,
  n_sim = PRED_ESS_N_SIM,
  trait_sd = PRED_ESS_TRAIT_SD,
  error_sd = PRED_ESS_ERROR_SD,
  seed = PRED_ESS_SEED + 2000
)
benchmark_summary_mono <- monotonize_benchmark_widths(benchmark_summary)

cat("Re-mapping only affected groups against extended benchmark...\n")
remapped_ess <- estimate_prediction_metric_ess(
  target_summary = affected_target_summary,
  benchmark_summary = benchmark_summary_mono,
  use_monotone = TRUE
)

if ("Random_ID" %in% names(affected_target_summary) &&
    nrow(affected_target_summary) == nrow(remapped_ess)) {
  remapped_ess$Random_ID <- affected_target_summary$Random_ID
}

existing_result <- if (file.exists(rds_file)) readRDS(rds_file) else list()
existing_ess_all <- existing_result$ess_summary_all
if (is.null(existing_ess_all)) {
  existing_ess_all <- read.csv(
    file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_all.csv"),
    stringsAsFactors = FALSE
  )
}

existing_keys <- make_key(existing_ess_all)
kept_ess <- existing_ess_all[!(existing_keys %in% affected_keys), , drop = FALSE]
ess_summary_all <- rbind(kept_ess, remapped_ess)

ess_summary_observed <- ess_summary_all[ess_summary_all$Subset_Type %in% c("dispersed", "clustered"), , drop = FALSE]
ess_summary_random <- ess_summary_all[ess_summary_all$Subset_Type == "random", , drop = FALSE]
ess_summary_vs_random <- summarize_piess_observed_vs_random(ess_summary_all)

ord <- c("N", "s", "Subset_Type", "Covariance_Model", "Covariance_Param", "Metric")
ess_summary_all <- ess_summary_all[do.call(order, ess_summary_all[intersect(ord, names(ess_summary_all))]), ]
ess_summary_observed <- ess_summary_observed[do.call(order, ess_summary_observed[intersect(ord, names(ess_summary_observed))]), ]
ess_summary_random <- ess_summary_random[do.call(order, ess_summary_random[intersect(ord, names(ess_summary_random))]), ]
ess_summary_vs_random <- ess_summary_vs_random[do.call(order, ess_summary_vs_random[intersect(ord, names(ess_summary_vs_random))]), ]

existing_result$benchmark_n <- benchmark_n_values
existing_result$benchmark_summary <- benchmark_summary
existing_result$benchmark_summary_mono <- benchmark_summary_mono
existing_result$target_summary <- target_summary
existing_result$ess_summary_all <- ess_summary_all
existing_result$ess_summary_observed <- ess_summary_observed
existing_result$ess_summary_random <- ess_summary_random
existing_result$ess_summary_vs_random <- ess_summary_vs_random

cat("Saving targeted remap results...\n")
saveRDS(existing_result, rds_file)
write.csv(benchmark_summary_mono, file.path(out_dir, "prediction_metric_ess_sensitivity_benchmark_summary.csv"), row.names = FALSE)
write.csv(ess_summary_all, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_all.csv"), row.names = FALSE)
write.csv(ess_summary_observed, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_observed.csv"), row.names = FALSE)
write.csv(ess_summary_random, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_random.csv"), row.names = FALSE)
write.csv(ess_summary_vs_random, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary_vs_random.csv"), row.names = FALSE)
write.csv(ess_summary_observed, file.path(out_dir, "prediction_metric_ess_sensitivity_ess_summary.csv"), row.names = FALSE)

cat("Done. Targeted extended-benchmark remap saved to:", out_dir, "\n")
