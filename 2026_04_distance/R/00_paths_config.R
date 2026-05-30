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

# Clustered subset exchange refinement parameters
CLUSTERED_MAX_EXCHANGE_ITERATIONS <- 100
CLUSTERED_EXCHANGE_TOL <- 1e-10
