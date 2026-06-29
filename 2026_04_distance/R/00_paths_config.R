# 00_paths_config.R
# Centralized path management and global parameters
# All scripts should source this file first.
# No setwd() calls anywhere in the project.

PROJECT_ROOT <- normalizePath(getwd())

RAW_DIR <- file.path(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
LOG_DIR <- file.path(PROJECT_ROOT, "logs")

TREE_FILE <- file.path(
  RAW_DIR,
  "MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"
)

TAXONOMY_FILE <- file.path(
  RAW_DIR,
  "taxonomy_mamPhy_5911species_toPublish.csv"
)

GLOBAL_SEED <- 20260428
N_NULL_REPS_DEFAULT <- 1000

CRICETIDAE_POOL_SIZES <- c(32, 64, 128, 256, 512)
SENSITIVITY_SUBSET_SIZES <- c(8, 16, 32, 64)

COV_SENSITIVITY_GRID <- data.frame(
  N = c(128, 128, 512, 512),
  s = c(8, 64, 8, 64)
)

LAMBDA_VALUES <- c(1.00, 0.75, 0.50, 0.25, 0.10, 0)
OU_HALF_LIFE_FRACS <- c(0.05, 0.10, 0.25, 0.50, 1.00)
EB_RATE_VALUES <- c(-4.0, -2.0, -1.0, -0.5, 0.0)

# Clustered subset exchange refinement parameters
CLUSTERED_MAX_EXCHANGE_ITERATIONS <- 10
CLUSTERED_EXCHANGE_TOL <- 1e-10

# ============================================================
# Prediction-metric-based effective sample size simulation
# ============================================================

PRED_ESS_POOL_NAME <- "C512"
PRED_ESS_N <- 512
PRED_ESS_S <- 64

# Main target subsets for prediction-metric ESS
# Main text should focus on dispersed and clustered.
PRED_ESS_TARGET_TYPES <- c("dispersed", "clustered")

# Optional random subsets may be summarized separately, but are not
# required for the main prediction-metric ESS table.
PRED_ESS_INCLUDE_RANDOM <- FALSE
PRED_ESS_RANDOM_N_SUMMARY <- 100

# Number of random subsets used for PIESS empirical significance.
# Use 100 for development / smoke tests; use 1000 for final manuscript
# if the table notes say "1000 random subsets".
PRED_ESS_RANDOM_N_FINAL <- 1000
PRED_ESS_RANDOM_N_DEV <- 100

# Default used by long-running PIESS scripts unless overridden by CLI.
PRED_ESS_RANDOM_N_DEFAULT <- PRED_ESS_RANDOM_N_SUMMARY

# Simulation settings
PRED_ESS_N_SIM <- 10000
PRED_ESS_BENCHMARK_N <- 4:32

# Two-scenario comparison only
PRED_ESS_LAMBDA_INDEPENDENT <- 0
PRED_ESS_LAMBDA_PHYLO <- 1

# Trait and prediction-error scale
PRED_ESS_TRAIT_SD <- 1
PRED_ESS_ERROR_VAR <- 0.1
PRED_ESS_ERROR_SD <- sqrt(PRED_ESS_ERROR_VAR)

# Numerical settings
PRED_ESS_SEED <- GLOBAL_SEED + 900000
PRED_ESS_MONOTONE_BENCHMARK <- TRUE
PRED_ESS_EIGEN_TOL <- 1e-10
