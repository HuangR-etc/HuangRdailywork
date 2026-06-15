# 04_run_cricetidae_sensitivity_grid.R
# Cricetidae sensitivity analysis: nested pools with varying N and s
# Nested pools: C32, C64, C128, C256, C512
# Subset sizes: s = 8, 16, 32, 64 (s < N)
# Each N-s combination runs independently
#
# Usage: Rscript scripts/04_run_cricetidae_sensitivity_grid.R
#        Rscript scripts/04_run_cricetidae_sensitivity_grid.R --N 128 --s 64 --n_null 100
args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_file[grepl("^--file=", args_file)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}
source("R/01_load_modules.R")
load_project_modules()

# ---- Parse optional command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)

RUN_N <- NA
RUN_S <- NA
N_NULL_REPS <- N_NULL_REPS_DEFAULT
OVERWRITE <- FALSE

i <- 1
while (i <= length(args)) {
  if (args[i] == "--N" && i < length(args)) {
    RUN_N <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--s" && i < length(args)) {
    RUN_S <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--n_null" && i < length(args)) {
    N_NULL_REPS <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--overwrite") {
    OVERWRITE <- TRUE
    i <- i + 1
  } else {
    i <- i + 1
  }
}

cat("=== 04_run_cricetidae_sensitivity_grid ===\n")
cat("N_NULL_REPS =", N_NULL_REPS, "\n")

# ---- Build grid ----
grid <- expand.grid(
  N = CRICETIDAE_POOL_SIZES,
  s = SENSITIVITY_SUBSET_SIZES
)
grid <- grid[grid$s < grid$N, ]
grid <- grid[order(grid$N, grid$s), ]

if (!is.na(RUN_N) && !is.na(RUN_S)) {
  grid <- grid[grid$N == RUN_N & grid$s == RUN_S, ]
  if (nrow(grid) == 0) {
    stop("No grid cell matches N = ", RUN_N, ", s = ", RUN_S)
  }
  cat("Running single case: N =", RUN_N, ", s =", RUN_S, "\n")
} else {
  cat("Running full grid:", nrow(grid), "cases\n")
}

# ---- Load nested pools ----
cat("Loading Cricetidae nested pools...\n")
nested_pools <- readRDS(file.path(PROCESSED_DIR, "cricetidae_nested_pools.rds"))

raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
sel_dir <- file.path(RESULTS_DIR, "sensitivity", "selected_subsets")
sum_dir <- file.path(RESULTS_DIR, "sensitivity", "summaries")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Run each case ----
all_dist_summaries <- list()
all_dep_summaries <- list()
all_selected_species <- list()

for (row_idx in seq_len(nrow(grid))) {
  N_i <- grid$N[row_idx]
  s_i <- grid$s[row_idx]
  
  pool_name <- paste0("C", N_i)
  pool_label <- paste0("Cricetidae_", pool_name)
  case_seed <- GLOBAL_SEED + N_i * 1000 + s_i
  
  cat("  Running: N =", N_i, ", s =", s_i,
      ", seed =", case_seed, "...\n")
  
  pool_tree <- nested_pools$pool_trees[[pool_name]]
  
  result <- run_empirical_case(
    pool_tree = pool_tree,
    pool_label = pool_label,
    subset_size = s_i,
    n_null_reps = N_NULL_REPS,
    seed = case_seed,
    out_dir = raw_dir,
    save_raw = TRUE,
    overwrite = OVERWRITE
  )
  
  all_dist_summaries[[length(all_dist_summaries) + 1]] <- result$distance_summary
  all_dep_summaries[[length(all_dep_summaries) + 1]] <- result$dependence_summary_BM
  all_selected_species[[length(all_selected_species) + 1]] <- result$selected_species
}

# ---- Save combined summaries ----
cat("Saving combined summaries...\n")

dist_combined <- do.call(rbind, all_dist_summaries)
write.csv(dist_combined,
          file.path(sum_dir, "cricetidae_sensitivity_distance_summary.csv"),
          row.names = FALSE)

dep_combined <- do.call(rbind, all_dep_summaries)
write.csv(dep_combined,
          file.path(sum_dir, "cricetidae_sensitivity_dependence_BM_summary.csv"),
          row.names = FALSE)

species_combined <- do.call(rbind, all_selected_species)
write.csv(species_combined,
          file.path(sel_dir, "cricetidae_sensitivity_selected_species.csv"),
          row.names = FALSE)

cat("Done. Cricetidae sensitivity results saved.\n")
