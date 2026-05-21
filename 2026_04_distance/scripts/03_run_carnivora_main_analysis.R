# 03_run_carnivora_main_analysis.R
# Carnivora empirical main analysis
# 261 species, s = 20, dispersed / clustered / random baseline
# Distance metrics + BM dependence diagnostics
#
# Usage: Rscript scripts/03_run_carnivora_main_analysis.R
setwd("/home/huangr/projects/2026_04_distance")
source("R/01_load_modules.R")
load_project_modules()

cat("=== 03_run_carnivora_main_analysis ===\n")

S <- 20
N_NULL <- 1000

out_dir <- file.path(RESULTS_DIR, "carnivora")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load Carnivora tree ----
cat("Loading Carnivora tree...\n")
carnivora_tree <- readRDS(file.path(PROCESSED_DIR, "carnivora_tree.rds"))
cat("  N =", ape::Ntip(carnivora_tree), "tips\n")

# ---- Run empirical case ----
cat("Running Carnivora main analysis (s =", S, ")...\n")
result <- run_empirical_case(
  pool_tree = carnivora_tree,
  pool_label = "Carnivora_N261",
  subset_size = S,
  n_null_reps = N_NULL,
  seed = GLOBAL_SEED,
  out_dir = out_dir,
  save_raw = TRUE,
  overwrite = FALSE
)

# ---- Save summary CSVs ----
cat("Saving summary CSVs...\n")

write.csv(result$distance_summary,
          file.path(out_dir, "carnivora_distance_summary_s20.csv"),
          row.names = FALSE)

write.csv(result$dependence_summary_BM,
          file.path(out_dir, "carnivora_dependence_BM_summary_s20.csv"),
          row.names = FALSE)

write.csv(result$selected_species,
          file.path(out_dir, "carnivora_selected_species_s20.csv"),
          row.names = FALSE)

# ---- Save raw random data ----
write.csv(result$random_dist,
          file.path(out_dir, "carnivora_random_distance_raw_s20.csv"),
          row.names = FALSE)

write.csv(result$random_dep_BM,
          file.path(out_dir, "carnivora_random_dependence_BM_raw_s20.csv"),
          row.names = FALSE)

# ---- Tree plots ----
cat("Saving tree plots...\n")

save_tree_plot_pdf(carnivora_tree, result$dispersed$final_subset_names,
                   file.path(out_dir, "carnivora_tree_dispersed_subset_s20.pdf"),
                   main = "Carnivora - Dispersed Subset (s=20)",
                   highlight_col = "red", cex = 0.5)

save_tree_plot_pdf(carnivora_tree, result$clustered$final_subset_names,
                   file.path(out_dir, "carnivora_tree_clustered_subset_s20.pdf"),
                   main = "Carnivora - Clustered Subset (s=20)",
                   highlight_col = "blue", cex = 0.5)

cat("Done. Carnivora results saved to:", out_dir, "\n")
