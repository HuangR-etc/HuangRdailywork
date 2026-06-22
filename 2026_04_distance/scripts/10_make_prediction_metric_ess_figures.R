# 10_make_prediction_metric_ess_figures.R
# Generate figures for prediction-metric-based effective sample size
#
# Produces a three-panel figure (RMSE, MAE, R^2) showing:
# - Independent benchmark curve (lambda = 0, n = 4:32)
# - Horizontal lines for lambda = 1 dispersed and clustered interval widths
# - Vertical dashed lines for interpolated prediction-metric ESS
#
# Usage: Rscript scripts/10_make_prediction_metric_ess_figures.R
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}
source("R/01_load_modules.R")
load_project_modules()

cat("=== 10_make_prediction_metric_ess_figures ===\n")

out_dir <- file.path(RESULTS_DIR, "prediction_metric_ess")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res_file <- file.path(out_dir, "prediction_metric_ess_C512_s64.rds")

if (!file.exists(res_file)) {
  stop("Missing result RDS: ", res_file,
       "\nRun scripts/09_run_prediction_metric_ess.R first.")
}

cat("Reading results from:", res_file, "\n")
res <- readRDS(res_file)

bench <- res$benchmark_summary_mono
target <- res$target_summary
ess <- res$ess_summary

metrics <- c("RMSE", "MAE", "R2")

# ---- PDF output ----
pdf_file <- file.path(out_dir, "Figure_prediction_metric_ess_C512_s64.pdf")
pdf(pdf_file, width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 0, 0))

for (m in metrics) {
  b <- bench[bench$Metric == m, ]
  t <- target[target$Metric == m, ]
  e <- ess[ess$Metric == m, ]
  
  # Determine y-axis limits
  y_range <- range(
    b$Interval_Width_95_Monotone,
    t$Interval_Width_95,
    na.rm = TRUE
  )
  y_pad <- diff(y_range) * 0.1
  y_lim <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  
  plot(b$Benchmark_N,
       b$Interval_Width_95_Monotone,
       type = "l",
       col = "black",
       lwd = 2,
       xlab = "Independent benchmark sample size",
       ylab = "95% empirical interval width",
       main = m,
       ylim = y_lim,
       cex.lab = 1.1,
       cex.main = 1.2)
  
  # Add horizontal lines for target subsets
  line_colors <- c("dispersed" = "#E41A1C", "clustered" = "#377EB8")
  line_lty <- c("dispersed" = 2, "clustered" = 3)
  
  for (stype in c("dispersed", "clustered")) {
    t_i <- t[t$Subset_Type == stype, ]
    e_i <- e[e$Subset_Type == stype, ]
    
    # Horizontal line for target interval width
    abline(h = t_i$Interval_Width_95,
           lty = line_lty[stype],
           col = line_colors[stype],
           lwd = 1.5)
    
    # Vertical line for ESS (if interpolated)
    if (is.finite(e_i$Prediction_Metric_ESS)) {
      abline(v = e_i$Prediction_Metric_ESS,
             lty = line_lty[stype],
             col = line_colors[stype],
             lwd = 1.5)
    }
  }
  
  # Add legend
  legend("topright",
         legend = c(expression(lambda == 0 ~ "benchmark"),
                    expression(lambda == 1 ~ "dispersed"),
                    expression(lambda == 1 ~ "clustered")),
         lty = c(1, 2, 3),
         col = c("black", "#E41A1C", "#377EB8"),
         lwd = c(2, 1.5, 1.5),
         bty = "n",
         cex = 0.9)
}

dev.off()
cat("Saved PDF:", pdf_file, "\n")
file.copy(pdf_file,
          file.path(out_dir, "Figure3_C512_s64_prediction_metric_ess.pdf"),
          overwrite = TRUE)
cat("Saved PDF:", file.path(out_dir, "Figure3_C512_s64_prediction_metric_ess.pdf"), "\n")

# ---- PNG output ----
png_file <- file.path(out_dir, "Figure_prediction_metric_ess_C512_s64.png")
png(png_file, width = 10, height = 4, units = "in", res = 300)
par(mfrow = c(1, 3), mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 0, 0))

for (m in metrics) {
  b <- bench[bench$Metric == m, ]
  t <- target[target$Metric == m, ]
  e <- ess[ess$Metric == m, ]
  
  y_range <- range(
    b$Interval_Width_95_Monotone,
    t$Interval_Width_95,
    na.rm = TRUE
  )
  y_pad <- diff(y_range) * 0.1
  y_lim <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  
  plot(b$Benchmark_N,
       b$Interval_Width_95_Monotone,
       type = "l",
       col = "black",
       lwd = 2,
       xlab = "Independent benchmark sample size",
       ylab = "95% empirical interval width",
       main = m,
       ylim = y_lim,
       cex.lab = 1.1,
       cex.main = 1.2)
  
  line_colors <- c("dispersed" = "#E41A1C", "clustered" = "#377EB8")
  line_lty <- c("dispersed" = 2, "clustered" = 3)
  
  for (stype in c("dispersed", "clustered")) {
    t_i <- t[t$Subset_Type == stype, ]
    e_i <- e[e$Subset_Type == stype, ]
    
    abline(h = t_i$Interval_Width_95,
           lty = line_lty[stype],
           col = line_colors[stype],
           lwd = 1.5)
    
    if (is.finite(e_i$Prediction_Metric_ESS)) {
      abline(v = e_i$Prediction_Metric_ESS,
             lty = line_lty[stype],
             col = line_colors[stype],
             lwd = 1.5)
    }
  }
  
  legend("topright",
         legend = c(expression(lambda == 0 ~ "benchmark"),
                    expression(lambda == 1 ~ "dispersed"),
                    expression(lambda == 1 ~ "clustered")),
         lty = c(1, 2, 3),
         col = c("black", "#E41A1C", "#377EB8"),
         lwd = c(2, 1.5, 1.5),
         bty = "n",
         cex = 0.9)
}

dev.off()
cat("Saved PNG:", png_file, "\n")
file.copy(png_file,
          file.path(out_dir, "Figure3_C512_s64_prediction_metric_ess.png"),
          overwrite = TRUE)
cat("Saved PNG:", file.path(out_dir, "Figure3_C512_s64_prediction_metric_ess.png"), "\n")

# ---- Also save a formatted table ----
cat("Saving formatted ESS table...\n")

table_out <- ess[, c("Subset_Type", "Metric", "Target_Interval_Width_95",
                      "Prediction_Metric_ESS_Label", "Match_Status")]
names(table_out) <- c("Subset_Type", "Metric", "95%_Interval_Width",
                       "Prediction_Metric_ESS", "Match_Status")

write.csv(table_out,
          file.path(out_dir, "Table_prediction_metric_ess_C512_s64.csv"),
          row.names = FALSE)

cat("Saved table:", file.path(out_dir, "Table_prediction_metric_ess_C512_s64.csv"), "\n")

# ---- Save supplement tables ----
cat("Saving supplement tables...\n")

# Supplement: benchmark curve
bench_out <- bench[, c("Benchmark_N", "Metric", "Interval_Width_95",
                        "Interval_Width_95_Monotone", "SD", "Q025", "Q500", "Q975")]
write.csv(bench_out,
          file.path(out_dir, "Supplement_prediction_metric_benchmark_curve.csv"),
          row.names = FALSE)

# Supplement: raw uncertainty summary
target_out <- target[, c("Subset_Type", "Lambda", "Metric", "Mean", "SD",
                          "Q025", "Q500", "Q975", "Interval_Width_95", "N_Sim")]
write.csv(target_out,
          file.path(out_dir, "Supplement_prediction_metric_raw_uncertainty_summary.csv"),
          row.names = FALSE)

cat("Done. All figures and tables saved to:", out_dir, "\n")
