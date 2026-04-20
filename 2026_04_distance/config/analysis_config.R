# Configuration for phylogenetic dispersed subset analysis
# This file contains all configurable parameters for the analysis

cfg <- list(
  # Random seed for reproducibility
  seed = 1,
  
  # Tree generation parameters
  tree_reps = 1,           # Number of tree replicates (default 1)
  large_n = 4096,           # Number of species in large trees
  small_n = 16,            # Number of species in small tree
  
  # Subset size parameters
  subset_large = 64,       # Subset size for large trees (256 species)
  subset_small = 4,        # Subset size for small tree (32 species)
  
  # Null distribution parameters
  null_reps_large = 1000,  # Number of random subsets for large trees
  null_reps_small = 1000,  # Number of random subsets for small tree
  exhaustive_small = TRUE, # Use exhaustive enumeration for small tree if TRUE
  
  # Trait simulation parameters
  trait_reps = 1,          # Number of trait simulation replicates (now 1 per tree)
  random_subset_reps_for_trait = 1000, # Number of random subsets for trait analysis (now 1000)
  
  # OU process parameters (alpha values for weak, moderate, strong)
  ou_alpha = c(0.2, 1, 5),
  
  # Lambda model parameters (lambda values for Pagel's lambda model)
  bm_lambda = c(1.0, 0.75,0.5),
  
  # Which results to run
  run_result1 = TRUE,      # Compare with random null
  run_result2 = TRUE,      # Test design components
  run_result3 = TRUE,      # Trait-level consequences
  run_result4 = TRUE,      # Heuristic vs exact optimum
  
  # Output directories
  output_dir = "outputs",
  tables_dir = "outputs/tables",
  figures_dir = "outputs/figures",
  rds_dir = "outputs/rds",
  logs_dir = "logs",
  
  # File naming conventions
  file_prefix = "phylo_disp",
  
  # Performance settings
  parallel_cores = 4,      # Number of cores for parallel processing
  save_intermediate = TRUE # Save intermediate results
)

# Function to create output directories if they don't exist
create_dirs <- function() {
  dirs <- c(cfg$output_dir, cfg$tables_dir, cfg$figures_dir, 
            cfg$rds_dir, cfg$logs_dir)
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

# Initialize directories
create_dirs()

# Set seed for reproducibility
set.seed(1)

# Print configuration summary
cat("Configuration loaded:\n")
cat("  Seed:", cfg$seed, "\n")
cat("  Tree replicates:", cfg$tree_reps, "\n")
cat("  Large tree size:", cfg$large_n, "\n")
cat("  Small tree size:", cfg$small_n, "\n")
cat("  Subset sizes:", cfg$subset_large, "(large),", cfg$subset_small, "(small)\n")
cat("  Null replicates:", cfg$null_reps_large, "(large),", cfg$null_reps_small, "(small)\n")
cat("  Trait replicates:", cfg$trait_reps, "\n")
cat("  Output directory:", cfg$output_dir, "\n")
