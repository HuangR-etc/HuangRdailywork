# 01_prepare_empirical_trees.R
# Prepare empirical trees: read mammal tree, extract Carnivora and Cricetidae,
# build nested Cricetidae pools.
#
# Usage: Rscript scripts/01_prepare_empirical_trees.R
setwd("/home/huangr/projects/2026_04_distance")
source("R/01_load_modules.R")
load_project_modules()

cat("=== 01_prepare_empirical_trees ===\n")

# ---- Read taxonomy ----
cat("Reading taxonomy...\n")
tax <- read.csv(TAXONOMY_FILE, stringsAsFactors = FALSE)

# ---- Read and clean mammal tree ----
cat("Reading and cleaning mammal tree...\n")
mammal_tree <- read_clean_mammal_tree(TREE_FILE)
cat("  Cleaned tree has", ape::Ntip(mammal_tree), "tips\n")

saveRDS(mammal_tree, file.path(PROCESSED_DIR, "mammal_4098_clean.rds"))

# ---- Extract Carnivora ----
cat("Extracting Carnivora...\n")
carnivora_tree <- extract_taxonomic_pool_tree(
  tree = mammal_tree,
  tax = tax,
  rank_col = "ord",
  rank_value = "Carnivora"
)
cat("  Carnivora tree has", ape::Ntip(carnivora_tree), "tips\n")

saveRDS(carnivora_tree, file.path(PROCESSED_DIR, "carnivora_tree.rds"))

# ---- Extract Cricetidae ----
cat("Extracting Cricetidae...\n")
cricetidae_tree <- extract_taxonomic_pool_tree(
  tree = mammal_tree,
  tax = tax,
  rank_col = "fam",
  rank_value = "CRICETIDAE"
)
cat("  Cricetidae tree has", ape::Ntip(cricetidae_tree), "tips\n")

saveRDS(cricetidae_tree, file.path(PROCESSED_DIR, "cricetidae_full_tree.rds"))

# ---- Build nested Cricetidae pools ----
cat("Building nested Cricetidae pools...\n")
nested_pools <- build_nested_candidate_pools(
  pool_tree = cricetidae_tree,
  pool_sizes = CRICETIDAE_POOL_SIZES,
  seed = GLOBAL_SEED
)

saveRDS(nested_pools, file.path(PROCESSED_DIR, "cricetidae_nested_pools.rds"))

# Save candidate pool CSVs
sel_dir <- file.path(RESULTS_DIR, "sensitivity", "selected_subsets")
dir.create(sel_dir, recursive = TRUE, showWarnings = FALSE)

for (pool_name in names(nested_pools$pools)) {
  csv_path <- file.path(sel_dir, paste0("Cricetidae_", pool_name, "_candidate_pool.csv"))
  write.csv(data.frame(Species = nested_pools$pools[[pool_name]]),
            csv_path, row.names = FALSE)
  cat("  Saved:", csv_path, "\n")
}

cat("Done. All empirical trees prepared.\n")
