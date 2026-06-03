# scripts/11_make_revised_empirical_figures.R
#
# Generate revised manuscript/supplement figures after replacing Carnivora
# with the Cricetidae C512 / s64 empirical case.
#
# This script DOES NOT rerun analyses. It reads the already computed raw case:
#   results/sensitivity/raw_cases/Cricetidae_C512_s64.rds
# and produces:
#   1. Supplementary full-page tree: dispersed subset on C512 tree
#   2. Supplementary full-page tree: clustered subset on C512 tree
#   3. Revised Figure 2abc: distance-metric random baseline densities
#   4. Revised Figure 3abcd: BM dependence diagnostics + lambda MeanESS
#   5. Traceability tables for Figure 2 and Figure 3 panels
#
# For panels where the dispersed subset is beyond all random baselines,
# the density curve is clipped to the observed-line boundary so the
# kernel-smoothed tail does not visually cross the line.
#
# Usage from the project root:
#   Rscript scripts/11_make_revised_empirical_figures.R
#
# If running elsewhere, first set your working directory to the project root.

source("R/01_load_modules.R")
load_project_modules()

cat("=== 11_make_revised_empirical_figures ===\n")

# ------------------------------------------------------------
# 0. Paths and inputs
# ------------------------------------------------------------
case_file <- file.path(
  RESULTS_DIR,
  "sensitivity", "raw_cases",
  "Cricetidae_C512_s64.rds"
)

if (!file.exists(case_file)) {
  stop(
    "Missing raw case file: ", case_file,
    "\nRun: Rscript scripts/04_run_cricetidae_sensitivity_grid.R --N 512 --s 64 --overwrite"
  )
}

out_dir <- file.path(RESULTS_DIR, "figures", "revised_empirical_C512_s64")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading raw case: ", case_file, "\n", sep = "")
case <- readRDS(case_file)

pool_tree <- case$pool_tree
subset_size <- case$s
pool_size <- case$N

if (pool_size != 512 || subset_size != 64) {
  warning(
    "Expected N = 512 and s = 64, but raw case has N = ", pool_size,
    " and s = ", subset_size, ". Figures will still be generated."
  )
}

disp_subset_names <- case$dispersed$final_subset_names
clust_subset_names <- case$clustered$final_subset_names

# Raw random distributions saved by run_empirical_case()
random_dist_metrics <- case$random_dist
random_dep_metrics <- case$random_dep_BM

# Recalculate observed metrics from the tree to avoid stale derived objects.
dist_obj <- create_distance_object(pool_tree)
disp_subset_idx <- match(disp_subset_names, dist_obj$tip_labels)
clust_subset_idx <- match(clust_subset_names, dist_obj$tip_labels)

if (anyNA(disp_subset_idx) || anyNA(clust_subset_idx)) {
  stop("Some selected subset names were not found in pool_tree$tip.label.")
}

disp_metrics <- calc_subset_metrics_extended(dist_obj$dist_mat, disp_subset_idx)
clust_metrics <- calc_subset_metrics_extended(dist_obj$dist_mat, clust_subset_idx)

V_bm <- make_bm_covariance(pool_tree)
disp_dep <- calc_dependence_from_V(V_bm, disp_subset_names)
clust_dep <- calc_dependence_from_V(V_bm, clust_subset_names)

# Save the exact plotting input summary for traceability.
plot_input_summary <- rbind(
  data.frame(
    Figure_Set = "Figure2_distance",
    Subset_Type = "dispersed",
    Metric = c("MinPD", "MeanPD", "MeanNND", "MaxPD"),
    Observed = c(disp_metrics$MinPD, disp_metrics$MeanPD,
                 disp_metrics$MeanNND, disp_metrics$MaxPD),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Figure_Set = "Figure2_distance",
    Subset_Type = "clustered",
    Metric = c("MinPD", "MeanPD", "MeanNND", "MaxPD"),
    Observed = c(clust_metrics$MinPD, clust_metrics$MeanPD,
                 clust_metrics$MeanNND, clust_metrics$MaxPD),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Figure_Set = "Figure3_dependence",
    Subset_Type = "dispersed",
    Metric = c("MeanOffCor", "MaxOffCor", "MeanESS"),
    Observed = c(disp_dep$off_mean, disp_dep$rmax, disp_dep$neff_mean),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Figure_Set = "Figure3_dependence",
    Subset_Type = "clustered",
    Metric = c("MeanOffCor", "MaxOffCor", "MeanESS"),
    Observed = c(clust_dep$off_mean, clust_dep$rmax, clust_dep$neff_mean),
    stringsAsFactors = FALSE
  )
)

write.csv(
  plot_input_summary,
  file.path(out_dir, "Revised_empirical_C512_s64_plot_input_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 1. Common plotting helpers
# ------------------------------------------------------------
expand_range <- function(x, add_frac = 0.08) {
  x <- x[is.finite(x)]
  xr <- range(x, na.rm = TRUE)
  dx <- diff(xr)
  if (!is.finite(dx) || dx == 0) {
    pad <- ifelse(abs(xr[1]) > 0, abs(xr[1]) * add_frac, 0.1)
  } else {
    pad <- dx * add_frac
  }
  c(xr[1] - pad, xr[2] + pad)
}

set_panel_par <- function() {
  par(
    mar = c(3.8, 4.0, 1.8, 0.8),
    mgp = c(2.0, 0.65, 0),
    tcl = -0.25,
    cex = 0.95
  )
}

safe_density <- function(x) {
  x <- x[is.finite(x)]
  if (length(unique(x)) < 2) return(NULL)
  density(x, na.rm = TRUE)
}

observed_range_status <- function(obs, null_values, favorable_direction = c("high", "low")) {
  favorable_direction <- match.arg(favorable_direction)
  null_values <- null_values[is.finite(null_values)]
  nmin <- min(null_values, na.rm = TRUE)
  nmax <- max(null_values, na.rm = TRUE)

  range_status <- if (!is.finite(obs)) {
    "NA"
  } else if (obs < nmin) {
    "below_random_min"
  } else if (obs > nmax) {
    "above_random_max"
  } else {
    "within_random_range"
  }

  more_extreme <- if (!is.finite(obs)) {
    NA
  } else if (favorable_direction == "high") {
    obs > nmax
  } else {
    obs < nmin
  }

  list(
    range_status = range_status,
    more_extreme_than_all_random = more_extreme,
    null_min = nmin,
    null_max = nmax
  )
}

make_one_panel_summary <- function(figure_set,
                                   panel,
                                   metric,
                                   null_values,
                                   obs_disp,
                                   obs_clust,
                                   favorable_direction_disp,
                                   favorable_direction_clust) {
  null_values <- null_values[is.finite(null_values)]
  disp_status <- observed_range_status(obs_disp, null_values, favorable_direction_disp)
  clust_status <- observed_range_status(obs_clust, null_values, favorable_direction_clust)

  # Directional empirical p-values, using the same +1 correction as the analysis.
  disp_p <- if (favorable_direction_disp == "high") calc_p_high(obs_disp, null_values) else calc_p_low(obs_disp, null_values)
  clust_p <- if (favorable_direction_clust == "high") calc_p_high(obs_clust, null_values) else calc_p_low(obs_clust, null_values)

  data.frame(
    Figure_Set = figure_set,
    Panel = panel,
    Metric = metric,
    Baseline_N = length(null_values),
    Baseline_Min = min(null_values, na.rm = TRUE),
    Baseline_Q025 = as.numeric(quantile(null_values, 0.025, na.rm = TRUE)),
    Baseline_Q25 = as.numeric(quantile(null_values, 0.25, na.rm = TRUE)),
    Baseline_Mean = mean(null_values, na.rm = TRUE),
    Baseline_Median = median(null_values, na.rm = TRUE),
    Baseline_Q75 = as.numeric(quantile(null_values, 0.75, na.rm = TRUE)),
    Baseline_Q975 = as.numeric(quantile(null_values, 0.975, na.rm = TRUE)),
    Baseline_Max = max(null_values, na.rm = TRUE),
    Baseline_SD = sd(null_values, na.rm = TRUE),
    Baseline_IQR = IQR(null_values, na.rm = TRUE),
    Dispersed_Observed = obs_disp,
    Dispersed_SES = calc_ses(obs_disp, null_values),
    Dispersed_P_value = disp_p,
    Dispersed_Favorable_Direction = favorable_direction_disp,
    Dispersed_Range_Status = disp_status$range_status,
    Dispersed_More_Extreme_Than_All_Random = disp_status$more_extreme_than_all_random,
    Clustered_Observed = obs_clust,
    Clustered_SES = calc_ses(obs_clust, null_values),
    Clustered_P_value = clust_p,
    Clustered_Favorable_Direction = favorable_direction_clust,
    Clustered_Range_Status = clust_status$range_status,
    Clustered_More_Extreme_Than_All_Random = clust_status$more_extreme_than_all_random,
    stringsAsFactors = FALSE
  )
}

figure2_3_summary <- do.call(rbind, list(
  make_one_panel_summary(
    figure_set = "Figure2_distance", panel = "2a", metric = "MinPD",
    null_values = random_dist_metrics$MinPD,
    obs_disp = disp_metrics$MinPD, obs_clust = clust_metrics$MinPD,
    favorable_direction_disp = "high", favorable_direction_clust = "low"
  ),
  make_one_panel_summary(
    figure_set = "Figure2_distance", panel = "2b", metric = "MeanPD",
    null_values = random_dist_metrics$MeanPD,
    obs_disp = disp_metrics$MeanPD, obs_clust = clust_metrics$MeanPD,
    favorable_direction_disp = "high", favorable_direction_clust = "low"
  ),
  make_one_panel_summary(
    figure_set = "Figure2_distance", panel = "2c", metric = "MeanNND",
    null_values = random_dist_metrics$MeanNND,
    obs_disp = disp_metrics$MeanNND, obs_clust = clust_metrics$MeanNND,
    favorable_direction_disp = "high", favorable_direction_clust = "low"
  ),
  make_one_panel_summary(
    figure_set = "Figure3_dependence", panel = "3a", metric = "MeanOffCor",
    null_values = random_dep_metrics$off_mean,
    obs_disp = disp_dep$off_mean, obs_clust = clust_dep$off_mean,
    favorable_direction_disp = "low", favorable_direction_clust = "high"
  ),
  make_one_panel_summary(
    figure_set = "Figure3_dependence", panel = "3b", metric = "MaxOffCor",
    null_values = random_dep_metrics$rmax,
    obs_disp = disp_dep$rmax, obs_clust = clust_dep$rmax,
    favorable_direction_disp = "low", favorable_direction_clust = "high"
  ),
  make_one_panel_summary(
    figure_set = "Figure3_dependence", panel = "3c", metric = "MeanESS",
    null_values = random_dep_metrics$neff_mean,
    obs_disp = disp_dep$neff_mean, obs_clust = clust_dep$neff_mean,
    favorable_direction_disp = "high", favorable_direction_clust = "low"
  )
))

write.csv(
  figure2_3_summary,
  file.path(out_dir, "Figure2_Figure3_baseline_observed_summary_C512_s64.csv"),
  row.names = FALSE
)

make_raw_plot_values <- function(figure_set, panel, metric, null_values, obs_disp, obs_clust) {
  null_df <- data.frame(
    Figure_Set = figure_set,
    Panel = panel,
    Metric = metric,
    Source = "random_baseline",
    Replicate_ID = seq_along(null_values),
    Value = as.numeric(null_values),
    stringsAsFactors = FALSE
  )
  obs_df <- data.frame(
    Figure_Set = figure_set,
    Panel = panel,
    Metric = metric,
    Source = c("dispersed_observed", "clustered_observed"),
    Replicate_ID = NA_integer_,
    Value = c(obs_disp, obs_clust),
    stringsAsFactors = FALSE
  )
  rbind(null_df, obs_df)
}

figure2_3_raw_values <- do.call(rbind, list(
  make_raw_plot_values("Figure2_distance", "2a", "MinPD", random_dist_metrics$MinPD, disp_metrics$MinPD, clust_metrics$MinPD),
  make_raw_plot_values("Figure2_distance", "2b", "MeanPD", random_dist_metrics$MeanPD, disp_metrics$MeanPD, clust_metrics$MeanPD),
  make_raw_plot_values("Figure2_distance", "2c", "MeanNND", random_dist_metrics$MeanNND, disp_metrics$MeanNND, clust_metrics$MeanNND),
  make_raw_plot_values("Figure3_dependence", "3a", "MeanOffCor", random_dep_metrics$off_mean, disp_dep$off_mean, clust_dep$off_mean),
  make_raw_plot_values("Figure3_dependence", "3b", "MaxOffCor", random_dep_metrics$rmax, disp_dep$rmax, clust_dep$rmax),
  make_raw_plot_values("Figure3_dependence", "3c", "MeanESS", random_dep_metrics$neff_mean, disp_dep$neff_mean, clust_dep$neff_mean)
))

write.csv(
  figure2_3_raw_values,
  file.path(out_dir, "Figure2_Figure3_plot_raw_values_C512_s64.csv"),
  row.names = FALSE
)

cat("\nRange checks for the panels you asked about:\n")
fig2b_row <- figure2_3_summary[figure2_3_summary$Panel == "2b", ]
fig3a_row <- figure2_3_summary[figure2_3_summary$Panel == "3a", ]
cat("  Figure 2b MeanPD: dispersed observed = ", fig2b_row$Dispersed_Observed,
    "; random max = ", fig2b_row$Baseline_Max,
    "; dispersed above all random = ", fig2b_row$Dispersed_More_Extreme_Than_All_Random, "\n", sep = "")
cat("  Figure 3a MeanOffCor: dispersed observed = ", fig3a_row$Dispersed_Observed,
    "; random min = ", fig3a_row$Baseline_Min,
    "; dispersed below all random = ", fig3a_row$Dispersed_More_Extreme_Than_All_Random, "\n\n", sep = "")

# Plot a random-baseline density panel.
# If the dispersed value is more extreme than every random baseline value,
# the x axis is locked at the dispersed value on the corresponding side.
# This prevents the Gaussian kernel tail from visually crossing the observed line.
plot_null_density_panel <- function(null_values,
                                    obs_disp = NULL,
                                    obs_clust = NULL,
                                    xlab,
                                    panel_label,
                                    show_legend = TRUE,
                                    legend_labels = NULL,
                                    favorable_direction_disp = c("high", "low"),
                                    lock_axis_if_disp_extreme = TRUE,
                                    hard_xlim = NULL) {
  favorable_direction_disp <- match.arg(favorable_direction_disp)
  null_values <- null_values[is.finite(null_values)]
  x_all <- c(null_values, obs_disp, obs_clust)

  if (!is.null(hard_xlim)) {
    xlim_use <- hard_xlim
  } else {
    xlim_use <- expand_range(x_all, add_frac = 0.10)

    if (lock_axis_if_disp_extreme && !is.null(obs_disp) && is.finite(obs_disp)) {
      null_min <- min(null_values, na.rm = TRUE)
      null_max <- max(null_values, na.rm = TRUE)
      if (favorable_direction_disp == "high" && obs_disp > null_max) {
        # Lock the right edge at the observed dispersed value.
        xlim_use[2] <- obs_disp
      }
      if (favorable_direction_disp == "low" && obs_disp < null_min) {
        # Lock the left edge at the observed dispersed value.
        xlim_use[1] <- obs_disp
      }
    }
  }

  den <- safe_density(null_values)

  if (is.null(den)) {
    # Fallback for near-constant null distributions.
    hist(
      null_values,
      breaks = 30,
      freq = FALSE,
      xlim = xlim_use,
      main = "",
      xlab = xlab,
      ylab = "Density",
      col = "grey90",
      border = "grey40"
    )
    ylim_top <- par("usr")[4]
  } else {
    # Clip the smoothed density curve to the locked x-axis range.
    # This avoids displaying kernel-smoothed probability mass beyond the observed line.
    keep <- den$x >= xlim_use[1] & den$x <= xlim_use[2]
    den_x <- den$x[keep]
    den_y <- den$y[keep]
    ylim_use <- c(0, max(den_y, na.rm = TRUE) * 1.10)

    plot(
      den_x,
      den_y,
      type = "l",
      xlim = xlim_use,
      ylim = ylim_use,
      main = "",
      xlab = xlab,
      ylab = "Density",
      lwd = 1.2,
      axes = FALSE
    )
    axis(1, cex.axis = 0.9)
    axis(2, las = 1, cex.axis = 0.9)
    box()
    ylim_top <- ylim_use[2]
  }

  # Solid = dispersed; dashed = clustered.
  if (!is.null(obs_disp) && is.finite(obs_disp)) {
    segments(obs_disp, 0, obs_disp, ylim_top, lwd = 1.8, lty = 1)
  }
  if (!is.null(obs_clust) && is.finite(obs_clust)) {
    segments(obs_clust, 0, obs_clust, ylim_top, lwd = 1.8, lty = 2)
  }

  mtext(panel_label, side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)

  if (show_legend) {
    if (is.null(legend_labels)) {
      legend_labels <- if (is.null(obs_clust)) "Dispersed" else c("Dispersed", "Clustered")
    }
    legend_lty <- if (length(legend_labels) == 1) 1 else c(1, 2)
    legend(
      "top",
      legend = legend_labels,
      lty = legend_lty,
      lwd = 2,
      bty = "n",
      cex = 0.85,
      seg.len = 4
    )
  }
}

# ------------------------------------------------------------
# 2. Full-page tree plots for supplementary figures
# ------------------------------------------------------------
plot_subset_on_full_tree <- function(tree,
                                     subset_names,
                                     file_name,
                                     main_title,
                                     point_cex = 0.75,
                                     outward_offset_frac = 0.007,
                                     page_width = 8.27,
                                     page_height = 11.69) {
  pdf(file_name, width = page_width, height = page_height, pointsize = 9)
  on.exit(dev.off(), add = TRUE)

  par(mar = c(0.2, 0.2, 1.0, 0.2), xpd = NA)

  plot.phylo(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    no.margin = TRUE,
    edge.color = "grey35",
    edge.width = 0.35
  )

  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tip <- ape::Ntip(tree)

  tip_idx <- match(subset_names, tree$tip.label)
  tip_idx <- tip_idx[is.finite(tip_idx)]

  tip_x <- pp$xx[tip_idx]
  tip_y <- pp$yy[tip_idx]

  root_x <- pp$xx[n_tip + 1]
  root_y <- pp$yy[n_tip + 1]

  dx <- tip_x - root_x
  dy <- tip_y - root_y
  r <- sqrt(dx^2 + dy^2)
  r[r == 0] <- 1

  all_tip_x <- pp$xx[seq_len(n_tip)]
  all_tip_y <- pp$yy[seq_len(n_tip)]
  all_r <- sqrt((all_tip_x - root_x)^2 + (all_tip_y - root_y)^2)
  offset <- max(all_r, na.rm = TRUE) * outward_offset_frac

  # Move selected-tip markers slightly outward so dots do not hide branch tips.
  tip_x_out <- tip_x + offset * dx / r
  tip_y_out <- tip_y + offset * dy / r

  points(
    tip_x_out,
    tip_y_out,
    pch = 21,
    bg = "black",
    col = "black",
    cex = point_cex,
    lwd = 0.35
  )

  mtext(main_title, side = 3, line = 0.2, font = 2, cex = 1.1)
}

plot_subset_on_full_tree(
  tree = pool_tree,
  subset_names = disp_subset_names,
  file_name = file.path(out_dir, "Supplement_Fig_Cricetidae_C512_s64_dispersed_tree_fullpage.pdf"),
  main_title = paste0("Cricetidae C512: phylogenetically dispersed subset (s = ", subset_size, ")")
)

plot_subset_on_full_tree(
  tree = pool_tree,
  subset_names = clust_subset_names,
  file_name = file.path(out_dir, "Supplement_Fig_Cricetidae_C512_s64_clustered_tree_fullpage.pdf"),
  main_title = paste0("Cricetidae C512: phylogenetically clustered subset (s = ", subset_size, ")")
)

# ------------------------------------------------------------
# 3. Revised Figure 2abc: distance metrics
#    Old Figure 2cde becomes new Figure 2abc.
# ------------------------------------------------------------
draw_fig2a_panel <- function() {
  plot_null_density_panel(
    null_values = random_dist_metrics$MinPD,
    obs_disp = disp_metrics$MinPD,
    obs_clust = NULL,
    xlab = "MinPD",
    panel_label = "a",
    show_legend = TRUE,
    legend_labels = "Dispersed",
    favorable_direction_disp = "high"
  )
}

draw_fig2b_panel <- function() {
  plot_null_density_panel(
    null_values = random_dist_metrics$MeanPD,
    obs_disp = disp_metrics$MeanPD,
    obs_clust = clust_metrics$MeanPD,
    xlab = "MeanPD",
    panel_label = "b",
    show_legend = TRUE,
    favorable_direction_disp = "high"
  )
}

draw_fig2c_panel <- function() {
  plot_null_density_panel(
    null_values = random_dist_metrics$MeanNND,
    obs_disp = disp_metrics$MeanNND,
    obs_clust = clust_metrics$MeanNND,
    xlab = "MeanNND",
    panel_label = "c",
    show_legend = TRUE,
    favorable_direction_disp = "high"
  )
}

pdf(
  file.path(out_dir, "Figure2abc_Cricetidae_C512_s64_distance_densities.pdf"),
  width = 12,
  height = 4.4,
  pointsize = 9
)
par(mfrow = c(1, 3))
set_panel_par()
draw_fig2a_panel()
draw_fig2b_panel()
draw_fig2c_panel()
dev.off()

# Optional single-panel exports, useful for editing in Illustrator/PowerPoint.
pdf(file.path(out_dir, "Figure2a_MinPD_density_C512_s64.pdf"), width = 3.8, height = 4.4, pointsize = 9)
set_panel_par(); draw_fig2a_panel(); dev.off()

pdf(file.path(out_dir, "Figure2b_MeanPD_density_C512_s64.pdf"), width = 3.8, height = 4.4, pointsize = 9)
set_panel_par(); draw_fig2b_panel(); dev.off()

pdf(file.path(out_dir, "Figure2c_MeanNND_density_C512_s64.pdf"), width = 3.8, height = 4.4, pointsize = 9)
set_panel_par(); draw_fig2c_panel(); dev.off()

# ------------------------------------------------------------
# 4. Revised Figure 3abcd: dependence diagnostics + lambda MeanESS
# ------------------------------------------------------------
# Use current project definition of MeanESS: 1' R^{-1} 1.
lambda_grid <- seq(0, 1, by = 0.1)

lambda_meanESS <- sapply(lambda_grid, function(lam) {
  V_lam <- lambda_transform_cov(V_bm, lam)
  calc_dependence_from_V(V_lam, disp_subset_names)$neff_mean
})

lambda_meanESS_df <- data.frame(
  lambda = lambda_grid,
  MeanESS = as.numeric(lambda_meanESS)
)

write.csv(
  lambda_meanESS_df,
  file.path(out_dir, "Figure3d_lambda_MeanESS_dispersed_C512_s64.csv"),
  row.names = FALSE
)

draw_fig3a_panel <- function() {
  plot_null_density_panel(
    null_values = random_dep_metrics$off_mean,
    obs_disp = disp_dep$off_mean,
    obs_clust = clust_dep$off_mean,
    xlab = "MeanOffCor",
    panel_label = "a",
    show_legend = TRUE,
    favorable_direction_disp = "low"
  )
}

draw_fig3b_panel <- function() {
  plot_null_density_panel(
    null_values = random_dep_metrics$rmax,
    obs_disp = disp_dep$rmax,
    obs_clust = clust_dep$rmax,
    xlab = "MaxOffCor",
    panel_label = "b",
    show_legend = TRUE,
    favorable_direction_disp = "low"
  )
}

draw_fig3c_panel <- function() {
  plot_null_density_panel(
    null_values = random_dep_metrics$neff_mean,
    obs_disp = disp_dep$neff_mean,
    obs_clust = clust_dep$neff_mean,
    xlab = "MeanESS",
    panel_label = "c",
    show_legend = TRUE,
    favorable_direction_disp = "high"
  )
}

draw_fig3d_panel <- function() {
  ylim_d <- expand_range(c(lambda_meanESS_df$MeanESS, subset_size), add_frac = 0.08)

  plot(
    lambda_meanESS_df$lambda,
    lambda_meanESS_df$MeanESS,
    type = "b",
    pch = 16,
    lwd = 1.4,
    xlim = c(0, 1),
    ylim = ylim_d,
    axes = FALSE,
    xlab = expression(lambda),
    ylab = "MeanESS",
    main = ""
  )

  axis(
    1,
    at = seq(0, 1, by = 0.1),
    labels = sprintf("%.1f", seq(0, 1, by = 0.1)),
    cex.axis = 0.9
  )
  axis(2, las = 1, cex.axis = 0.9)
  box()

  # Optional reference line for nominal subset size.
  # Uncomment if you want to show the nominal s = 64 ceiling visually.
  # abline(h = subset_size, lty = 3, lwd = 1.2)

  mtext("d", side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
}

pdf(
  file.path(out_dir, "Figure3abcd_Cricetidae_C512_s64_dependence_densities_lambda.pdf"),
  width = 14.4,
  height = 4.4,
  pointsize = 9
)
par(mfrow = c(1, 4))
set_panel_par()
draw_fig3a_panel()
draw_fig3b_panel()
draw_fig3c_panel()
draw_fig3d_panel()
dev.off()

# Optional single-panel exports.
pdf(file.path(out_dir, "Figure3a_MeanOffCor_density_C512_s64.pdf"), width = 3.2, height = 4.4, pointsize = 9)
set_panel_par(); draw_fig3a_panel(); dev.off()

pdf(file.path(out_dir, "Figure3b_MaxOffCor_density_C512_s64.pdf"), width = 3.2, height = 4.4, pointsize = 9)
set_panel_par(); draw_fig3b_panel(); dev.off()

pdf(file.path(out_dir, "Figure3c_MeanESS_density_C512_s64.pdf"), width = 3.2, height = 4.4, pointsize = 9)
set_panel_par(); draw_fig3c_panel(); dev.off()

pdf(file.path(out_dir, "Figure3d_lambda_MeanESS_C512_s64.pdf"), width = 3.2, height = 4.4, pointsize = 9)
set_panel_par(); draw_fig3d_panel(); dev.off()

# ------------------------------------------------------------
# 5. Final console report
# ------------------------------------------------------------
cat("\nDone. Revised empirical figure PDFs saved to:\n")
cat(out_dir, "\n\n")
cat("Main outputs:\n")
cat("  Supplement_Fig_Cricetidae_C512_s64_dispersed_tree_fullpage.pdf\n")
cat("  Supplement_Fig_Cricetidae_C512_s64_clustered_tree_fullpage.pdf\n")
cat("  Figure2abc_Cricetidae_C512_s64_distance_densities.pdf\n")
cat("  Figure3abcd_Cricetidae_C512_s64_dependence_densities_lambda.pdf\n")
cat("\nTraceability tables:\n")
cat("  Revised_empirical_C512_s64_plot_input_summary.csv\n")
cat("  Figure2_Figure3_baseline_observed_summary_C512_s64.csv\n")
cat("  Figure2_Figure3_plot_raw_values_C512_s64.csv\n")
cat("  Figure3d_lambda_MeanESS_dispersed_C512_s64.csv\n")
