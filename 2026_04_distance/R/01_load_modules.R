# 01_load_modules.R
# Unified package loading and project module sourcing

load_project_modules <- function() {
  suppressPackageStartupMessages({
    library(ape)
  })
  
  source("R/00_paths_config.R")
  source("R/02_tree_io.R")
  source("R/03_distance_metrics.R")
  source("R/04_subset_dispersed.R")
  source("R/05_subset_clustered.R")
  source("R/06_random_baseline.R")
  source("R/07_dependence_diagnostics.R")
  source("R/08_covariance_models.R")
  source("R/09_summary_tables.R")
  source("R/10_plotting_trees.R")
  source("R/11_sanity_tree_helpers.R")
  source("R/12_prediction_metric_ess.R")
}
