# Example: Running Result 1 Random-Tree Analysis
# This script demonstrates how to use the new random-tree module

# Set working directory to project root
setwd("/home/huangr/projects/2026_04_distance")

# Load required libraries
library(ape)
library(phytools)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load existing modules (needed for the analysis functions)
source("R/distance_metrics.R")
source("R/objective_compare.R")
source("R/subset_greedy.R")
source("R/subset_exchange.R")
source("R/subset_random.R")

# Load the new random-tree module
source("R/result1_random_trees.R")

cat(paste0("\n", strrep("=", 80), "\n"))
cat("EXAMPLE: RESULT 1 RANDOM-TREE ANALYSIS\n")
cat(paste0(strrep("=", 80), "\n"))

# Option 1: Run a quick test with small parameters
cat("\nOption 1: Running test analysis (small parameters for speed)...\n")
cat("This will generate 5 random trees with 50 tips each.\n")
cat("Each tree will have subset size 10 with 100 null replicates.\n\n")

test_results <- test_result1_random_trees()

cat("\nTest analysis completed.\n")
cat("Output directory: test_outputs/result1_random_trees_test\n")

# Option 2: Run full analysis with default parameters
cat(paste0("\n", strrep("-", 80), "\n"))
cat("\nOption 2: Full analysis with default parameters...\n")
cat("Default parameters:\n")
cat("  Number of random trees: 100\n")
cat("  Tips per tree: 256\n")
cat("  Subset size: 20\n")
cat("  Null replicates per tree: 1000\n")
cat("  Global seed: 20260409\n")
cat("  Output directory: outputs/result1_random_trees\n\n")

cat("To run the full analysis, uncomment the following lines:\n")
cat("```r\n")
cat("# full_results <- run_result1_random_trees()\n")
cat("```\n")

# Option 3: Run with custom parameters
cat(paste0("\n", strrep("-", 80), "\n"))
cat("\nOption 3: Custom analysis...\n")
cat("Example with custom parameters:\n")
cat("```r\n")
cat("custom_results <- run_result1_random_trees(\n")
cat("  n_trees = 50,           # Fewer trees for faster analysis\n")
cat("  n_tips = 128,           # Smaller trees\n")
cat("  subset_size = 15,       # Different subset size\n")
cat("  n_null_reps = 500,      # Fewer null replicates\n")
cat("  global_seed = 12345,    # Custom seed\n")
cat("  output_dir = \"outputs/result1_random_trees_custom\"\n")
cat(")\n")
cat("```\n")

# Option 4: Analyze individual components
cat(paste0("\n", strrep("-", 80), "\n"))
cat("\nOption 4: Using individual functions...\n")

# Generate random trees
cat("\n1. Generating random trees:\n")
cat("```r\n")
cat("trees <- generate_random_trees(n_trees = 10, n_tips = 100)\n")
cat("```\n")

# Analyze a single tree
cat("\n2. Analyzing a single tree:\n")
cat("```r\n")
cat("tree_result <- analyze_single_random_tree(trees[[1]], \n")
cat("  subset_size = 20, n_null_reps = 100)\n")
cat("```\n")

# Create directory structure
cat("\n3. Creating output directories:\n")
cat("```r\n")
cat("dirs <- create_random_tree_dirs(\"my_outputs/random_trees\")\n")
cat("```\n")

# Create cross-tree summary from existing data
cat("\n4. Creating cross-tree summary:\n")
cat("```r\n")
cat("# Assuming you have a per_tree_summary.csv file\n")
cat("per_tree_df <- read.csv(\"outputs/result1_random_trees/summary/per_tree_summary.csv\")\n")
cat("cross_summary <- create_cross_tree_summary(per_tree_df)\n")
cat("```\n")

# Visualize results
cat("\n5. Creating visualizations:\n")
cat("```r\n")
cat("# Create visualizations from existing results\n")
cat("create_random_tree_visualizations(per_tree_df, all_results, \"my_figures\")\n")
cat("```\n")

# Show example of accessing results
cat(paste0("\n", strrep("-", 80), "\n"))
cat("\nAccessing Results:\n")

cat("\nAfter running the analysis, you can access:\n")
cat("1. Per-tree summary: results$per_tree_summary\n")
cat("2. Cross-tree summary: results$cross_tree_summary\n")
cat("3. Configuration: results$config\n")
cat("4. Individual tree results: results$all_results[[1]]\n")

cat("\nExample:\n")
cat("```r\n")
cat("# Get median z-score for MinPD\n")
cat("median_z_minpd <- median(results$per_tree_summary$z_MinPD)\n")
cat("cat(\"Median z-score for MinPD:\", median_z_minpd, \"\\n\")\n")
cat("\n")
cat("# Get proportion of trees where observed exceeds 95th percentile\n")
cat("prop_above95 <- mean(results$per_tree_summary$above95_MinPD)\n")
cat("cat(\"Proportion above 95th percentile:\", prop_above95, \"\\n\")\n")
cat("```\n")

# Key outputs
cat(paste0("\n", strrep("-", 80), "\n"))
cat("\nKey Outputs:\n")

cat("\nThe analysis generates:\n")
cat("1. Per-tree results in outputs/result1_random_trees/per_tree/\n")
cat("   - Each tree has its own directory with tree.rds, observed_subset.csv, etc.\n")
cat("\n2. Summary files in outputs/result1_random_trees/summary/\n")
cat("   - per_tree_summary.csv: One row per tree with all metrics\n")
cat("   - cross_tree_summary.csv: Summary statistics across all trees\n")
cat("   - random_tree_config.json: Configuration used\n")
cat("\n3. Figures in outputs/result1_random_trees/figures/\n")
cat("   - zscore_boxplot.pdf: Distribution of z-scores across trees\n")
cat("   - percentile_boxplot.pdf: Distribution of percentiles\n")
cat("   - above95_barplot.pdf: Proportion of trees above 95th percentile\n")
cat("   - example_tree_null_dist.pdf: Example tree null distribution\n")
cat("   - compare_with_existing_result1.pdf: Comparison with fixed trees\n")
cat("\n4. Report in outputs/result1_random_trees/summary/random_tree_analysis_report.txt\n")

cat(paste0("\n", strrep("=", 80), "\n"))
cat("SUMMARY\n")
cat(paste0(strrep("=", 80), "\n"))

cat("\nThe Result 1 Random-Tree module provides:\n")
cat("1. Independent analysis that doesn't interfere with existing Result 1\n")
cat("2. Robustness testing across random phylogenetic topologies\n")
cat("3. Standardized metrics (z-scores, percentiles) for cross-tree comparison\n")
cat("4. Comprehensive visualizations and summary statistics\n")
cat("5. Comparison with original fixed tree results\n")

cat("\nThis module addresses the limitation of fixed tree examples by showing\n")
cat("that the method's effectiveness generalizes to random tree topologies.\n")

cat("\nFor more details, see:\n")
cat("- R/result1_random_trees.R: Main implementation\n")
cat("- test_random_trees.R: Test script\n")
cat("- README.md: Updated documentation\n")
