# scripts/07_make_metric_tables.R
# Build manuscript-ready metric-specific tables from existing summary CSVs.
#
# Usage:
#   Rscript scripts/07_make_metric_tables.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}

source("R/01_load_modules.R")
load_project_modules()

cat("=== 07_make_metric_tables ===\n")

library(stats)

table_dir <- file.path(RESULTS_DIR, "figures", "metric_tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read input summaries ----
sens_dist <- read.csv(file.path(
  RESULTS_DIR, "sensitivity", "summaries",
  "cricetidae_sensitivity_distance_summary.csv"
), stringsAsFactors = FALSE)

sens_dep <- read.csv(file.path(
  RESULTS_DIR, "sensitivity", "summaries",
  "cricetidae_sensitivity_dependence_BM_summary.csv"
), stringsAsFactors = FALSE)

cov_sens <- read.csv(file.path(
  RESULTS_DIR, "covariance_sensitivity",
  "cricetidae_covariance_sensitivity_summary.csv"
), stringsAsFactors = FALSE)

# ---- Constants ----
N_values <- c(32, 64, 128, 256, 512)
s_values <- c(8, 16, 32, 64)
subset_types <- c("dispersed", "clustered")

# ============================================================
# 3. Helper formatting functions
# ============================================================

# 3.1 p-value formatting
format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# 3.2 Numeric value formatting
format_value <- function(x, digits = 2) {
  if (is.na(x)) return("NA")
  sprintf(paste0("%.", digits, "f"), x)
}

# 3.3 Value + p formatting
format_value_p <- function(value, p, digits = 2) {
  paste0(format_value(value, digits), " (", format_p(p), ")")
}

# 3.4 Significance stars
p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

# ============================================================
# 4. N–s metric table function
# ============================================================
# Used for: distance summary, BM dependence summary
# Produces format like screenshot 1.

make_ns_metric_table <- function(df,
                                 metric_name,
                                 value_digits = 2,
                                 N_values = c(32, 64, 128, 256, 512),
                                 s_values = c(8, 16, 32, 64),
                                 subset_types = c("dispersed", "clustered")) {
  
  out <- expand.grid(
    Subset_Type = subset_types,
    s = s_values,
    stringsAsFactors = FALSE
  )
  
  # Keep desired order
  out$Subset_Type <- factor(out$Subset_Type, levels = subset_types)
  out <- out[order(out$Subset_Type, out$s), ]
  out$Subset_Type <- as.character(out$Subset_Type)
  
  for (N_i in N_values) {
    col_name <- paste0("N_", N_i)
    out[[col_name]] <- "—"
    
    for (row_i in seq_len(nrow(out))) {
      stype_i <- out$Subset_Type[row_i]
      s_i <- out$s[row_i]
      
      if (s_i >= N_i) {
        out[[col_name]][row_i] <- "—"
        next
      }
      
      hit <- df[
        df$Metric == metric_name &
          df$Subset_Type == stype_i &
          df$N == N_i &
          df$s == s_i,
      ]
      
      if (nrow(hit) == 0) {
        out[[col_name]][row_i] <- "—"
      } else {
        out[[col_name]][row_i] <- format_value_p(
          value = hit$Observed[1],
          p = hit$P_value[1],
          digits = value_digits
        )
      }
    }
  }
  
  names(out)[names(out) == "Subset_Type"] <- "Subset type"
  out
}

# ============================================================
# 5. Batch generate distance and BM-dependence tables
# ============================================================

# 5.1 Distance metric tables
distance_metrics <- c("MinPD", "MeanPD", "MeanNND", "MaxPD")

for (metric_i in distance_metrics) {
  cat("Building distance table:", metric_i, "\n")
  
  tab <- make_ns_metric_table(
    df = sens_dist,
    metric_name = metric_i,
    value_digits = 2
  )
  
  write.csv(
    tab,
    file.path(table_dir, paste0("Table_distance_", metric_i, ".csv")),
    row.names = FALSE
  )
}

# 5.2 BM dependence metric tables
dependence_metrics <- c("off_mean", "rmax", "neff_mean")

for (metric_i in dependence_metrics) {
  cat("Building BM dependence table:", metric_i, "\n")
  
  digits_i <- ifelse(metric_i == "neff_mean", 2, 3)
  
  tab <- make_ns_metric_table(
    df = sens_dep,
    metric_name = metric_i,
    value_digits = digits_i
  )
  
  write.csv(
    tab,
    file.path(table_dir, paste0("Table_BM_dependence_", metric_i, ".csv")),
    row.names = FALSE
  )
}

# ============================================================
# 5b. Manuscript-style N–s metric tables: observed value + stars
# ============================================================
# These process tables use the same N–s layout as the manuscript tables.
#
# Cell format:
#   observed_value*
#   observed_value**
#   observed_value***
#
# No exact p-values are shown in the cells.
# Non-significant cells show the observed value without stars.
# Dashes indicate combinations not analyzed because s >= N.

# ------------------------------------------------------------
# 5b.1 Value + significance-star formatting
# ------------------------------------------------------------

format_value_stars <- function(value, p, digits = 2) {
  paste0(
    format_value(value, digits),
    p_to_stars(p)
  )
}

# ------------------------------------------------------------
# 5b.2 N–s metric table function: value + stars
# ------------------------------------------------------------

make_ns_metric_value_star_table <- function(df,
                                            metric_name,
                                            value_digits = 2,
                                            N_values = c(32, 64, 128, 256, 512),
                                            s_values = c(8, 16, 32, 64),
                                            subset_types = c("dispersed", "clustered")) {
  
  out <- expand.grid(
    Subset_Type = subset_types,
    s = s_values,
    stringsAsFactors = FALSE
  )
  
  # Keep same row order as the existing N-s tables
  out$Subset_Type <- factor(out$Subset_Type, levels = subset_types)
  out <- out[order(out$Subset_Type, out$s), ]
  out$Subset_Type <- as.character(out$Subset_Type)
  
  for (N_i in N_values) {
    col_name <- paste0("N_", N_i)
    out[[col_name]] <- "—"
    
    for (row_i in seq_len(nrow(out))) {
      stype_i <- out$Subset_Type[row_i]
      s_i <- out$s[row_i]
      
      if (s_i >= N_i) {
        out[[col_name]][row_i] <- "—"
        next
      }
      
      hit <- df[
        df$Metric == metric_name &
          df$Subset_Type == stype_i &
          df$N == N_i &
          df$s == s_i,
      ]
      
      if (nrow(hit) == 0) {
        out[[col_name]][row_i] <- "—"
      } else {
        out[[col_name]][row_i] <- format_value_stars(
          value = hit$Observed[1],
          p = hit$P_value[1],
          digits = value_digits
        )
      }
    }
  }
  
  names(out)[names(out) == "Subset_Type"] <- "Subset type"
  names(out)[names(out) == "s"] <- "Subset size"
  
  out
}

# ------------------------------------------------------------
# 5b.3 Output distance metric tables in value + stars format
# ------------------------------------------------------------
# Distance metrics are used as same-format supplementary/process tables.

for (metric_i in distance_metrics) {
  cat("Building value-star N-s distance table:", metric_i, "\n")
  
  tab_star <- make_ns_metric_value_star_table(
    df = sens_dist,
    metric_name = metric_i,
    value_digits = 2,
    N_values = N_values,
    s_values = s_values,
    subset_types = subset_types
  )
  
  write.csv(
    tab_star,
    file.path(
      table_dir,
      paste0("Table_distance_", metric_i, "_value_stars.csv")
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 5b.4 Output BM dependence metric tables in value + stars format
# ------------------------------------------------------------
# off_mean, rmax, and neff_mean correspond to same-format
# supplementary BM-dependence tables.
#
# Use 3 decimals for correlations and 2 decimals for MeanESS.

for (metric_i in dependence_metrics) {
  cat("Building value-star N-s BM dependence table:", metric_i, "\n")
  
  digits_i <- ifelse(metric_i == "neff_mean", 2, 3)
  
  tab_star <- make_ns_metric_value_star_table(
    df = sens_dep,
    metric_name = metric_i,
    value_digits = digits_i,
    N_values = N_values,
    s_values = s_values,
    subset_types = subset_types
  )
  
  write.csv(
    tab_star,
    file.path(
      table_dir,
      paste0("Table_BM_dependence_", metric_i, "_value_stars.csv")
    ),
    row.names = FALSE
  )
}
# ============================================================
# 6. Covariance panel A / B table function (value + p)
# ============================================================
# Panel A: lambda_BM
# Panel B: OU
# Each cell: D: value (p) / C: value (p)

make_covariance_panel_table <- function(cov_df,
                                        metric_name,
                                        model_name,
                                        value_digits = 2) {
  
  df <- cov_df[
    cov_df$Metric == metric_name &
      cov_df$Covariance_Model == model_name,
  ]
  
  if (nrow(df) == 0) {
    stop("No rows found for metric = ", metric_name,
         " and model = ", model_name)
  }
  
  # Define covariance settings
  param_values <- sort(unique(df$Covariance_Param))
  
  if (model_name == "lambda_BM") {
    setting_labels <- paste0("\u03bb = ", sprintf("%.2f", param_values))
  } else if (model_name == "OU") {
    setting_labels <- paste0("h / H = ", sprintf("%.2f", param_values))
  } else {
    setting_labels <- as.character(param_values)
  }
  
  out <- data.frame(
    Covariance_setting = setting_labels,
    stringsAsFactors = FALSE
  )
  
  case_grid <- unique(df[, c("N", "s")])
  case_grid <- case_grid[order(case_grid$N, case_grid$s), ]
  
  for (case_i in seq_len(nrow(case_grid))) {
    N_i <- case_grid$N[case_i]
    s_i <- case_grid$s[case_i]
    
    col_name <- paste0("N_", N_i, "_s_", s_i)
    out[[col_name]] <- NA_character_
    
    for (p_i in seq_along(param_values)) {
      param_i <- param_values[p_i]
      
      d_row <- df[
        df$N == N_i &
          df$s == s_i &
          df$Covariance_Param == param_i &
          df$Subset_Type == "dispersed",
      ]
      
      c_row <- df[
        df$N == N_i &
          df$s == s_i &
          df$Covariance_Param == param_i &
          df$Subset_Type == "clustered",
      ]
      
      if (nrow(d_row) == 0 || nrow(c_row) == 0) {
        out[[col_name]][p_i] <- "—"
      } else {
        d_txt <- paste0(
          "D: ",
          format_value(d_row$Observed[1], value_digits),
          " (",
          format_p(d_row$P_value[1]),
          ")"
        )
        
        c_txt <- paste0(
          "C: ",
          format_value(c_row$Observed[1], value_digits),
          " (",
          format_p(c_row$P_value[1]),
          ")"
        )
        
        out[[col_name]][p_i] <- paste0(d_txt, " / ", c_txt)
      }
    }
  }
  
  names(out)[1] <- "Covariance setting"
  out
}

# ============================================================
# 7. Covariance value-star compact table function
# ============================================================
# Produces manuscript format like:
#   7.82** / 4.13***
# where the left value is dispersed and the right value is clustered.
# No explicit D/C labels and no exact p-values are shown.

make_covariance_panel_value_star_table <- function(cov_df,
                                                   metric_name,
                                                   model_name,
                                                   value_digits = 2) {
  
  df <- cov_df[
    cov_df$Metric == metric_name &
      cov_df$Covariance_Model == model_name,
  ]
  
  if (nrow(df) == 0) {
    stop("No rows found for metric = ", metric_name,
         " and model = ", model_name)
  }
  
  # Ensure numeric ordering of covariance settings
  param_values <- sort(as.numeric(unique(df$Covariance_Param)))
  
  if (model_name == "lambda_BM") {
    setting_labels <- paste0("\u03bb = ", sprintf("%.2f", param_values))
  } else if (model_name == "OU") {
    setting_labels <- paste0("h / H = ", sprintf("%.2f", param_values))
  } else {
    setting_labels <- as.character(param_values)
  }
  
  out <- data.frame(
    Covariance_setting = setting_labels,
    stringsAsFactors = FALSE
  )
  
  case_grid <- unique(df[, c("N", "s")])
  case_grid <- case_grid[order(case_grid$N, case_grid$s), ]
  
  for (case_i in seq_len(nrow(case_grid))) {
    N_i <- case_grid$N[case_i]
    s_i <- case_grid$s[case_i]
    
    col_name <- paste0("N_", N_i, "_s_", s_i)
    out[[col_name]] <- NA_character_
    
    for (p_i in seq_along(param_values)) {
      param_i <- param_values[p_i]
      
      d_row <- df[
        df$N == N_i &
          df$s == s_i &
          as.numeric(df$Covariance_Param) == param_i &
          df$Subset_Type == "dispersed",
      ]
      
      c_row <- df[
        df$N == N_i &
          df$s == s_i &
          as.numeric(df$Covariance_Param) == param_i &
          df$Subset_Type == "clustered",
      ]
      
      if (nrow(d_row) == 0 || nrow(c_row) == 0) {
        out[[col_name]][p_i] <- "—"
      } else {
        d_txt <- paste0(
          format_value(d_row$Observed[1], value_digits),
          p_to_stars(d_row$P_value[1])
        )
        
        c_txt <- paste0(
          format_value(c_row$Observed[1], value_digits),
          p_to_stars(c_row$P_value[1])
        )
        
        out[[col_name]][p_i] <- paste0(d_txt, " / ", c_txt)
      }
    }
  }
  
  names(out)[1] <- "Covariance setting"
  out
}
# ============================================================
# 8. Batch output covariance panel tables
# ============================================================

for (metric_i in dependence_metrics) {
  cat("Building covariance panel tables:", metric_i, "\n")
  
  digits_i <- ifelse(metric_i == "neff_mean", 2, 3)
  
  # ----------------------------------------------------------
  # 8.1 Full value + p tables
  # These are useful for supplementary tables or checking.
  # Format: D: value (p) / C: value (p)
  # ----------------------------------------------------------
  
  tab_lambda_value_p <- make_covariance_panel_table(
    cov_df = cov_sens,
    metric_name = metric_i,
    model_name = "lambda_BM",
    value_digits = digits_i
  )
  
  write.csv(
    tab_lambda_value_p,
    file.path(
      table_dir,
      paste0("Table_covariance_", metric_i, "_PanelA_lambdaBM_value_p.csv")
    ),
    row.names = FALSE
  )
  
  tab_ou_value_p <- make_covariance_panel_table(
    cov_df = cov_sens,
    metric_name = metric_i,
    model_name = "OU",
    value_digits = digits_i
  )
  
  write.csv(
    tab_ou_value_p,
    file.path(
      table_dir,
      paste0("Table_covariance_", metric_i, "_PanelB_OU_value_p.csv")
    ),
    row.names = FALSE
  )
  
  tab_eb_value_p <- make_covariance_panel_table(
    cov_df = cov_sens,
    metric_name = metric_i,
    model_name = "EB",
    value_digits = digits_i
  )
  
  write.csv(
    tab_eb_value_p,
    file.path(
      table_dir,
      paste0("Table_covariance_", metric_i, "_PanelC_EB_value_p.csv")
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 8.2 Manuscript compact value + star tables
  # Format: dispersed_value* / clustered_value***
  # No D/C labels and no exact p-values.
  # ----------------------------------------------------------
  
  tab_lambda_value_stars <- make_covariance_panel_value_star_table(
    cov_df = cov_sens,
    metric_name = metric_i,
    model_name = "lambda_BM",
    value_digits = digits_i
  )
  
  write.csv(
    tab_lambda_value_stars,
    file.path(
      table_dir,
      paste0("Table_covariance_", metric_i, "_PanelA_lambdaBM_value_stars.csv")
    ),
    row.names = FALSE
  )
  
  tab_ou_value_stars <- make_covariance_panel_value_star_table(
    cov_df = cov_sens,
    metric_name = metric_i,
    model_name = "OU",
    value_digits = digits_i
  )
  
  write.csv(
    tab_ou_value_stars,
    file.path(
      table_dir,
      paste0("Table_covariance_", metric_i, "_PanelB_OU_value_stars.csv")
    ),
    row.names = FALSE
  )
  
  tab_eb_value_stars <- make_covariance_panel_value_star_table(
    cov_df = cov_sens,
    metric_name = metric_i,
    model_name = "EB",
    value_digits = digits_i
  )
  
  write.csv(
    tab_eb_value_stars,
    file.path(
      table_dir,
      paste0("Table_covariance_", metric_i, "_PanelC_EB_value_stars.csv")
    ),
    row.names = FALSE
  )

  if (metric_i == "neff_mean") {
    write.csv(
      tab_lambda_value_stars,
      file.path(table_dir, "Table3_MIESS_covariance_PanelA_lambdaBM_value_stars.csv"),
      row.names = FALSE
    )
    write.csv(
      tab_ou_value_stars,
      file.path(table_dir, "Table3_MIESS_covariance_PanelB_OU_value_stars.csv"),
      row.names = FALSE
    )
    write.csv(
      tab_eb_value_stars,
      file.path(table_dir, "Table3_MIESS_covariance_PanelC_EB_value_stars.csv"),
      row.names = FALSE
    )
  }
}
# ============================================================
# 8b. P-value-only transparency tables
# ============================================================
# These tables are intended as supplementary process data.
#
# For supplementary/process N–s tables:
#   Keep the same row/column structure as make_ns_metric_table().
#   Rows: subset type × subset size s
#   Columns: candidate-pool size N
#   Cells: empirical one-sided p-values only.
#
# For covariance-panel supplementary/process tables:
#   Because dispersed and clustered values are combined in one cell
#   in the manuscript-style tables, p-value-only versions are split
#   into separate dispersed and clustered tables.
#
# All cells contain raw p-values only:
#   no observed values, no stars, no D/C labels.

# ------------------------------------------------------------
# 8b.1 Raw p-value formatting
# ------------------------------------------------------------
# Use this instead of format_p(), because these transparency tables
# should show the original p-values rather than thresholds such as <0.001.

format_p_raw <- function(p) {
  if (length(p) == 0 || is.na(p)) return("NA")
  format(p, scientific = FALSE, trim = TRUE, digits = 15)
}

# ------------------------------------------------------------
# 8b.2 N–s p-value-only table function
# ------------------------------------------------------------
# Used for N–s supplementary/process tables.
# Keeps dispersed and clustered in the same table as separate row blocks.

make_ns_metric_pvalue_table <- function(df,
                                        metric_name,
                                        N_values = c(32, 64, 128, 256, 512),
                                        s_values = c(8, 16, 32, 64),
                                        subset_types = c("dispersed", "clustered")) {
  
  out <- expand.grid(
    Subset_Type = subset_types,
    s = s_values,
    stringsAsFactors = FALSE
  )
  
  # Keep the same row order as make_ns_metric_table()
  out$Subset_Type <- factor(out$Subset_Type, levels = subset_types)
  out <- out[order(out$Subset_Type, out$s), ]
  out$Subset_Type <- as.character(out$Subset_Type)
  
  for (N_i in N_values) {
    col_name <- paste0("N_", N_i)
    out[[col_name]] <- "—"
    
    for (row_i in seq_len(nrow(out))) {
      stype_i <- out$Subset_Type[row_i]
      s_i <- out$s[row_i]
      
      if (s_i >= N_i) {
        out[[col_name]][row_i] <- "—"
        next
      }
      
      hit <- df[
        df$Metric == metric_name &
          df$Subset_Type == stype_i &
          df$N == N_i &
          df$s == s_i,
      ]
      
      if (nrow(hit) == 0) {
        out[[col_name]][row_i] <- "—"
      } else {
        out[[col_name]][row_i] <- format_p_raw(hit$P_value[1])
      }
    }
  }
  
  names(out)[names(out) == "Subset_Type"] <- "Subset type"
  names(out)[names(out) == "s"] <- "Subset size"
  
  out
}

# ------------------------------------------------------------
# 8b.3 Covariance p-value-only panel table function
# ------------------------------------------------------------
# Used for covariance supplementary/process tables.
# Dispersed and clustered are generated as separate tables.

make_covariance_panel_pvalue_table <- function(cov_df,
                                               metric_name,
                                               model_name,
                                               subset_type) {
  
  df <- cov_df[
    cov_df$Metric == metric_name &
      cov_df$Covariance_Model == model_name &
      cov_df$Subset_Type == subset_type,
  ]
  
  if (nrow(df) == 0) {
    stop("No rows found for metric = ", metric_name,
         ", model = ", model_name,
         ", and subset type = ", subset_type)
  }
  
  param_values <- sort(as.numeric(unique(df$Covariance_Param)))
  
  if (model_name == "lambda_BM") {
    setting_labels <- paste0("\u03bb = ", sprintf("%.2f", param_values))
  } else if (model_name == "OU") {
    setting_labels <- paste0("h / H = ", sprintf("%.2f", param_values))
  } else {
    setting_labels <- as.character(param_values)
  }
  
  out <- data.frame(
    Covariance_setting = setting_labels,
    stringsAsFactors = FALSE
  )
  
  case_grid <- unique(df[, c("N", "s")])
  case_grid <- case_grid[order(case_grid$N, case_grid$s), ]
  
  for (case_i in seq_len(nrow(case_grid))) {
    N_i <- case_grid$N[case_i]
    s_i <- case_grid$s[case_i]
    
    col_name <- paste0("N_", N_i, "_s_", s_i)
    out[[col_name]] <- "—"
    
    for (p_i in seq_along(param_values)) {
      param_i <- param_values[p_i]
      
      hit <- df[
        df$N == N_i &
          df$s == s_i &
          as.numeric(df$Covariance_Param) == param_i,
      ]
      
      if (nrow(hit) == 0) {
        out[[col_name]][p_i] <- "—"
      } else {
        out[[col_name]][p_i] <- format_p_raw(hit$P_value[1])
      }
    }
  }
  
  names(out)[1] <- "Covariance setting"
  out
}

# ------------------------------------------------------------
# 8b.4 Output p-value-only tables for N–s supplementary/process tables
# ------------------------------------------------------------
# These use the same N–s layout as the manuscript/supplementary
# N–s metric tables.
#
# Distance metrics:
#   Distance metrics correspond to supplementary/process distance tables.

for (metric_i in distance_metrics) {
  cat("Building p-value-only N-s distance table:", metric_i, "\n")
  
  tab_p <- make_ns_metric_pvalue_table(
    df = sens_dist,
    metric_name = metric_i,
    N_values = N_values,
    s_values = s_values,
    subset_types = subset_types
  )
  
  write.csv(
    tab_p,
    file.path(
      table_dir,
      paste0("Pvalues_distance_", metric_i, ".csv")
    ),
    row.names = FALSE
  )
}

# BM dependence metrics:
#   off_mean, rmax, and neff_mean correspond to supplementary
#   BM dependence tables.

for (metric_i in dependence_metrics) {
  cat("Building p-value-only N-s BM dependence table:", metric_i, "\n")
  
  tab_p <- make_ns_metric_pvalue_table(
    df = sens_dep,
    metric_name = metric_i,
    N_values = N_values,
    s_values = s_values,
    subset_types = subset_types
  )
  
  write.csv(
    tab_p,
    file.path(
      table_dir,
      paste0("Pvalues_BM_dependence_", metric_i, ".csv")
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 8b.5 Output p-value-only tables for covariance supplementary/process tables
# ------------------------------------------------------------
# These use the same covariance-panel layout as the manuscript-style
# covariance tables, but dispersed and clustered are split.
#
# off_mean, rmax, and neff_mean correspond to covariance supplementary/process tables.

for (metric_i in dependence_metrics) {
  cat("Building p-value-only covariance panel tables:", metric_i, "\n")
  
  for (model_i in c("lambda_BM", "OU", "EB")) {
    
    panel_label <- ifelse(model_i == "lambda_BM", "PanelA_lambdaBM",
                          ifelse(model_i == "OU", "PanelB_OU", "PanelC_EB"))
    
    for (stype_i in subset_types) {
      
      tab_p <- make_covariance_panel_pvalue_table(
        cov_df = cov_sens,
        metric_name = metric_i,
        model_name = model_i,
        subset_type = stype_i
      )
      
      write.csv(
        tab_p,
        file.path(
          table_dir,
          paste0(
            "Pvalues_covariance_",
            metric_i, "_",
            panel_label, "_",
            stype_i,
            ".csv"
          )
        ),
        row.names = FALSE
      )
    }
  }
}

cat("P-value-only transparency tables saved to:\n")
cat(table_dir, "\n")

# ============================================================
# 10. Main Table 1 and Supplementary Tables S5-S6: BM dependence N-s tables
#     Format: observed value + significance stars
# ============================================================
# S5: MeanOffCor = off_mean
# S6: MaxOffCor  = rmax
# Table 1: MIESS = neff_mean
#
# Cell format:
#   observed_value*
#   observed_value**
#   observed_value***
#
# Non-significant cells show only the observed value.
# En dashes indicate combinations not analyzed because s >= N.

# ------------------------------------------------------------
# 10.1 Helper: value + stars
# ------------------------------------------------------------

format_value_stars <- function(value, p, digits = 2) {
  paste0(
    format_value(value, digits),
    p_to_stars(p)
  )
}

# ------------------------------------------------------------
# 10.2 Helper: make N-s dependence table in manuscript format
# ------------------------------------------------------------

make_bm_dependence_supp_table <- function(df,
                                          metric_name,
                                          value_digits = 3,
                                          N_values = c(32, 64, 128, 256, 512),
                                          s_values = c(8, 16, 32, 64),
                                          subset_types = c("dispersed", "clustered")) {
  
  out <- expand.grid(
    `Subset type` = subset_types,
    s = s_values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Keep row order: dispersed 8/16/32/64, then clustered 8/16/32/64
  out$`Subset type` <- factor(out$`Subset type`, levels = subset_types)
  out <- out[order(out$`Subset type`, out$s), ]
  out$`Subset type` <- as.character(out$`Subset type`)
  
  for (N_i in N_values) {
    col_name <- paste0("N = ", N_i)
    out[[col_name]] <- "—"
    
    for (row_i in seq_len(nrow(out))) {
      stype_i <- out$`Subset type`[row_i]
      s_i <- out$s[row_i]
      
      if (s_i >= N_i) {
        out[[col_name]][row_i] <- "—"
        next
      }
      
      hit <- df[
        df$Metric == metric_name &
          df$Subset_Type == stype_i &
          df$N == N_i &
          df$s == s_i,
      ]
      
      if (nrow(hit) == 0) {
        out[[col_name]][row_i] <- "—"
      } else {
        out[[col_name]][row_i] <- format_value_stars(
          value = hit$Observed[1],
          p = hit$P_value[1],
          digits = value_digits
        )
      }
    }
  }
  
  out
}

# ------------------------------------------------------------
# 10.3 Table metadata
# ------------------------------------------------------------

bm_supp_tables <- data.frame(
  Filename = c(
    "Supplementary_Table_S5_BM_MeanOffCor_value_stars.csv",
    "Supplementary_Table_S6_BM_MaxOffCor_value_stars.csv",
    "Table1_MIESS_BM_nested_value_stars.csv"
  ),
  Metric = c("off_mean", "rmax", "neff_mean"),
  Digits = c(3, 3, 2),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(bm_supp_tables))) {
  tab <- make_bm_dependence_supp_table(
    df = sens_dep,
    metric_name = bm_supp_tables$Metric[i],
    value_digits = bm_supp_tables$Digits[i],
    N_values = N_values,
    s_values = s_values,
    subset_types = subset_types
  )
  write.csv(tab, file.path(table_dir, bm_supp_tables$Filename[i]), row.names = FALSE)
}

cat("BM dependence manuscript tables saved to:\n")
cat(table_dir, "\n")
