# 13_make_prediction_metric_ess_sensitivity_tables.R
# Build manuscript-ready PIESS and direct predictive-metric shift tables.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
source("R/01_load_modules.R")
load_project_modules()

cat("=== 13_make_prediction_metric_ess_sensitivity_tables ===\n")

dash <- "\u2014"
out_dir <- file.path(RESULTS_DIR, "figures", "metric_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    warning("Missing input, skipping related tables: ", path)
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_table <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE, na = "")
  cat("  Saved:", filename, "\n")
}

make_piess_ns_value_star_table <- function(df, metric_name,
                                           N_values = c(32, 64, 128, 256, 512),
                                           s_values = c(8, 16, 32, 64),
                                           subset_types = c("dispersed", "clustered"),
                                           value_col = "PIESS_Display_Stars") {
  out <- expand.grid(Subset_type = subset_types, Subset_size = s_values, stringsAsFactors = FALSE)
  for (N_i in N_values) out[[paste0("N_", N_i)]] <- dash
  hit_all <- df[df$Metric == metric_name & df$Condition == "BM_lambda1_target", , drop = FALSE]
  for (i in seq_len(nrow(out))) {
    for (N_i in N_values) {
      col_i <- paste0("N_", N_i)
      if (out$Subset_size[i] >= N_i) next
      hit <- hit_all[hit_all$Subset_Type == out$Subset_type[i] &
                       hit_all$s == out$Subset_size[i] & hit_all$N == N_i, , drop = FALSE]
      if (nrow(hit) > 0 && value_col %in% names(hit)) {
        val <- hit[[value_col]][1]
        if (!is.na(val) && nzchar(as.character(val))) out[[col_i]][i] <- as.character(val)
      }
    }
  }
  out
}

panel_spec <- function(model_name) {
  if (model_name == "lambda_BM") {
    return(list(panel = "PanelA_lambdaBM",
                labels = sprintf("\u03bb = %.2f", c(0, 0.10, 0.25, 0.50, 0.75, 1.00)),
                params = c(0, 0.10, 0.25, 0.50, 0.75, 1.00)))
  }
  if (model_name == "OU") {
    return(list(panel = "PanelB_OU",
                labels = sprintf("h/H = %.2f", c(0.05, 0.10, 0.25, 0.50, 1.00)),
                params = c(0.05, 0.10, 0.25, 0.50, 1.00)))
  }
  if (model_name == "EB") {
    return(list(panel = "PanelC_EB",
                labels = sprintf("EB r = %.2f", c(-0.10, -0.50, -1.00)),
                params = c(-0.10, -0.50, -1.00)))
  }
  stop("Unknown model_name: ", model_name)
}

case_columns <- data.frame(
  N = c(128, 128, 512, 512),
  s = c(8, 64, 8, 64),
  col = c("N_128_s_8", "N_128_s_64", "N_512_s_8", "N_512_s_64"),
  stringsAsFactors = FALSE
)

make_piess_covariance_panel_value_star_table <- function(df, metric_name, model_name) {
  spec <- panel_spec(model_name)
  out <- data.frame("Covariance setting" = spec$labels, stringsAsFactors = FALSE, check.names = FALSE)
  for (col_i in case_columns$col) out[[col_i]] <- dash
  hit_all <- df[df$Metric == metric_name & df$Covariance_Model == model_name, , drop = FALSE]
  for (i in seq_len(nrow(out))) {
    for (j in seq_len(nrow(case_columns))) {
      vals <- c()
      for (stype in c("dispersed", "clustered")) {
        hit <- hit_all[hit_all$N == case_columns$N[j] & hit_all$s == case_columns$s[j] &
                         hit_all$Subset_Type == stype &
                         abs(hit_all$Covariance_Param - spec$params[i]) < 1e-8, , drop = FALSE]
        vals <- c(vals, if (nrow(hit) > 0) hit$PIESS_Display_Stars[1] else dash)
      }
      out[[case_columns$col[j]]][i] <- paste(vals, collapse = " / ")
    }
  }
  out
}

make_piess_covariance_pvalue_table <- function(df, metric_name, model_name, subset_type) {
  spec <- panel_spec(model_name)
  out <- data.frame("Covariance setting" = spec$labels, stringsAsFactors = FALSE, check.names = FALSE)
  for (col_i in case_columns$col) out[[col_i]] <- dash
  hit_all <- df[df$Metric == metric_name & df$Covariance_Model == model_name &
                  df$Subset_Type == subset_type, , drop = FALSE]
  for (i in seq_len(nrow(out))) {
    for (j in seq_len(nrow(case_columns))) {
      hit <- hit_all[hit_all$N == case_columns$N[j] & hit_all$s == case_columns$s[j] &
                       abs(hit_all$Covariance_Param - spec$params[i]) < 1e-8, , drop = FALSE]
      if (nrow(hit) > 0 && is.finite(hit$P_value[1])) {
        out[[case_columns$col[j]]][i] <- sprintf("%.6g", hit$P_value[1])
      }
    }
  }
  out
}

make_direct_shift_covariance_panel_table <- function(df, metric_name, model_name) {
  spec <- panel_spec(model_name)
  out <- data.frame("Covariance setting" = spec$labels, stringsAsFactors = FALSE, check.names = FALSE)
  for (col_i in case_columns$col) out[[col_i]] <- dash
  hit_all <- df[df$Metric == metric_name & df$Covariance_Model == model_name &
                  df$Subset_Type == "dispersed", , drop = FALSE]
  for (i in seq_len(nrow(out))) {
    for (j in seq_len(nrow(case_columns))) {
      hit <- hit_all[hit_all$N == case_columns$N[j] & hit_all$s == case_columns$s[j] &
                       abs(hit_all$Covariance_Param - spec$params[i]) < 1e-8, , drop = FALSE]
      if (nrow(hit) > 0 && is.finite(hit$Percent_Change[1])) {
        out[[case_columns$col[j]]][i] <- sprintf("%.2f", hit$Percent_Change[1])
      }
    }
  }
  out
}

full_grid_vs_random <- read_optional_csv(file.path(
  RESULTS_DIR, "prediction_metric_ess_full_grid", "prediction_metric_ess_full_grid_ess_summary_vs_random.csv"
))
full_grid_shift <- read_optional_csv(file.path(
  RESULTS_DIR, "prediction_metric_ess_full_grid", "prediction_metric_direct_shift_full_grid.csv"
))
cov_vs_random <- read_optional_csv(file.path(
  RESULTS_DIR, "prediction_metric_ess_sensitivity", "prediction_metric_ess_sensitivity_ess_summary_vs_random.csv"
))
cov_shift <- read_optional_csv(file.path(
  RESULTS_DIR, "prediction_metric_ess_sensitivity", "prediction_metric_direct_shift_covariance_sensitivity.csv"
))

if (!is.null(full_grid_vs_random)) {
  piess_files <- list(
    R2 = "Table3_PIESS_R2_BM_nested_value_stars.csv",
    RMSE = "TableS7_PIESS_RMSE_BM_nested_value_stars.csv",
    MAE = "TableS8_PIESS_MAE_BM_nested_value_stars.csv"
  )
  for (metric_i in names(piess_files)) {
    write_table(make_piess_ns_value_star_table(full_grid_vs_random, metric_i), piess_files[[metric_i]])
    write_table(
      make_piess_ns_value_star_table(full_grid_vs_random, metric_i, value_col = "P_value"),
      paste0("Pvalues_PIESS_", metric_i, "_BM_nested.csv")
    )
  }
}

if (!is.null(cov_vs_random)) {
  prefix_map <- list(
    R2 = "Table5_PIESS_R2_covariance",
    RMSE = "TableS11_PIESS_RMSE_covariance",
    MAE = "TableS12_PIESS_MAE_covariance"
  )
  for (metric_i in names(prefix_map)) {
    for (model_i in c("lambda_BM", "OU", "EB")) {
      spec <- panel_spec(model_i)
      write_table(
        make_piess_covariance_panel_value_star_table(cov_vs_random, metric_i, model_i),
        paste0(prefix_map[[metric_i]], "_", spec$panel, "_value_stars.csv")
      )
      for (stype in c("dispersed", "clustered")) {
        write_table(
          make_piess_covariance_pvalue_table(cov_vs_random, metric_i, model_i, stype),
          paste0("Pvalues_PIESS_", metric_i, "_covariance_", spec$panel, "_", stype, ".csv")
        )
      }
    }
  }
}

if (!is.null(full_grid_shift)) {
  tab1 <- full_grid_shift[
    full_grid_shift$N == 512 & full_grid_shift$s == 64 &
      full_grid_shift$Subset_Type == "dispersed" &
      full_grid_shift$Condition == "BM_lambda1_target" &
      full_grid_shift$Metric %in% c("RMSE", "MAE", "R2"),
    , drop = FALSE
  ]
  if (nrow(tab1) > 0) {
    tab1 <- tab1[match(c("RMSE", "MAE", "R2"), tab1$Metric), , drop = FALSE]
    tab1 <- data.frame(
      Metric = tab1$Metric,
      Independent_Mean = round(tab1$Independent_Mean, 4),
      BM_Lambda1_Mean = round(tab1$Scenario_Mean, 4),
      Absolute_Change = round(tab1$Absolute_Change, 4),
      Percent_Change = round(tab1$Percent_Change, 2),
      stringsAsFactors = FALSE
    )
    write_table(tab1, "Table1_direct_metric_shift_C512_s64_dispersed.csv")
  } else {
    warning("No rows found for Table1 direct metric shift.")
  }
}

if (!is.null(cov_shift)) {
  shift_prefix <- list(
    R2 = "Table6_R2_percent_inflation_covariance",
    RMSE = "TableS13_RMSE_percent_deflation_covariance",
    MAE = "TableS14_MAE_percent_deflation_covariance"
  )
  for (metric_i in names(shift_prefix)) {
    for (model_i in c("lambda_BM", "OU", "EB")) {
      spec <- panel_spec(model_i)
      write_table(
        make_direct_shift_covariance_panel_table(cov_shift, metric_i, model_i),
        paste0(shift_prefix[[metric_i]], "_", spec$panel, ".csv")
      )
    }
  }
}

cat("Done. Metric tables saved to:", out_dir, "\n")

# ---- Extra Table 1-style direct-shift tables for subset/error-variance variants ----
make_table1_direct_shift <- function(direct_shift, subset_type, filename) {
  tab <- direct_shift[
    direct_shift$N == 512 & direct_shift$s == 64 &
      direct_shift$Subset_Type == subset_type &
      direct_shift$Condition == "BM_lambda1_target" &
      direct_shift$Metric %in% c("RMSE", "MAE", "R2"),
    , drop = FALSE
  ]
  if (nrow(tab) == 0) {
    warning("No rows found for Table 1-style direct metric shift: ", filename)
    return(invisible(FALSE))
  }
  tab <- tab[match(c("RMSE", "MAE", "R2"), tab$Metric), , drop = FALSE]
  tab <- data.frame(
    Metric = tab$Metric,
    Independent_Mean = round(tab$Independent_Mean, 4),
    BM_Lambda1_Mean = round(tab$Scenario_Mean, 4),
    Absolute_Change = round(tab$Absolute_Change, 4),
    Percent_Change = round(tab$Percent_Change, 2),
    stringsAsFactors = FALSE
  )
  write_table(tab, filename)
  invisible(TRUE)
}

write_table1_variant_set <- function(input_dir, suffix = "") {
  f <- file.path(input_dir, "prediction_metric_direct_shift_full_grid.csv")
  if (!file.exists(f)) {
    warning("Missing direct-shift input for Table 1-style variants: ", f)
    return(invisible(FALSE))
  }
  direct_shift <- read.csv(f, stringsAsFactors = FALSE)
  make_table1_direct_shift(
    direct_shift,
    "dispersed",
    paste0("Table1_direct_metric_shift_C512_s64_dispersed", suffix, ".csv")
  )
  make_table1_direct_shift(
    direct_shift,
    "clustered",
    paste0("Table1_direct_metric_shift_C512_s64_clustered", suffix, ".csv")
  )
  invisible(TRUE)
}

write_table1_variant_set(
  file.path(RESULTS_DIR, "prediction_metric_ess_full_grid"),
  suffix = ""
)
write_table1_variant_set(
  file.path(RESULTS_DIR, "prediction_metric_ess_full_grid_error_var_0p1"),
  suffix = "_error_var_0p1"
)
