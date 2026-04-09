# Test script to verify the implementation
# This script runs basic tests on each module

cat(paste0("\n", strrep("=", 80), "\n"))
cat("TESTING PHYLOGENETIC DISPERSED SUBSET IMPLEMENTATION\n")
cat(paste0(strrep("=", 80), "\n"))

# Set working directory
setwd("/home/huangr/projects/2026_04_distance")

# Load configuration
cat("\n1. Loading configuration...\n")
source("config/analysis_config.R")

# Create a copy of cfg and modify for testing
test_cfg <- cfg
test_cfg$large_n <- 32  # Must be power of 2 for balanced tree
test_cfg$small_n <- 16  # Must be power of 2 for balanced tree
test_cfg$subset_large <- 5
test_cfg$subset_small <- 3
test_cfg$null_reps_large <- 50
test_cfg$null_reps_small <- 50
test_cfg$trait_reps <- 5
test_cfg$random_subset_reps_for_trait <- 10
test_cfg$exhaustive_small <- FALSE

cat("Test configuration:\n")
print(test_cfg)

# Set seed
set.seed(test_cfg$seed)

# Test 1: Tree generation
cat("\n2. Testing tree generation...\n")
source("R/tree_generators.R")

trees <- generate_all_trees(test_cfg)
cat("Generated trees:", names(trees), "\n")
cat("Large balanced tree tips:", length(trees$large_balanced$tip.label), "\n")
cat("Large ladder tree tips:", length(trees$large_ladder$tip.label), "\n")
cat("Small balanced tree tips:", length(trees$small_balanced$tip.label), "\n")

# Test 2: Distance metrics
cat("\n3. Testing distance metrics...\n")
source("R/distance_metrics.R")

# Create a distance object
dist_obj <- create_distance_object(trees$small_balanced)
cat("Distance matrix dimensions:", dim(dist_obj$dist_mat), "\n")
cat("Tip labels:", length(dist_obj$tip_labels), "\n")

# Test subset metrics
test_subset <- c(1, 5, 10, 12, 15)
metrics <- calc_subset_metrics(dist_obj$dist_mat, test_subset)
cat("Test subset metrics:\n")
print(metrics)

# Test 3: Objective comparison
cat("\n4. Testing objective comparison...\n")
source("R/objective_compare.R")

# Create two metric sets
metrics1 <- list(MinPD = 10, MeanPD = 20, MeanNND = 5)
metrics2 <- list(MinPD = 9, MeanPD = 25, MeanNND = 6)

cat("Metrics 1: MinPD=", metrics1$MinPD, "MeanPD=", metrics1$MeanPD, "MeanNND=", metrics1$MeanNND, "\n")
cat("Metrics 2: MinPD=", metrics2$MinPD, "MeanPD=", metrics2$MeanPD, "MeanNND=", metrics2$MeanNND, "\n")
cat("Is metrics1 better than metrics2?", is_better_lexico_max(metrics1, metrics2), "\n")
cat("Is metrics2 better than metrics1?", is_better_lexico_max(metrics2, metrics1), "\n")

# Test 4: Greedy algorithm
cat("\n5. Testing greedy algorithm...\n")
source("R/subset_greedy.R")

greedy_result <- build_subset_greedy(dist_obj, subset_size = 5, maximize = TRUE)
cat("Greedy algorithm completed\n")
cat("Selected subset size:", length(greedy_result$subset), "\n")
cat("Selected subset names:", greedy_result$subset_names, "\n")
cat("Final metrics:\n")
print(greedy_result$metrics)

# Test 5: Exchange refinement
cat("\n6. Testing exchange refinement...\n")
source("R/subset_exchange.R")

exchange_result <- refine_subset_exchange(dist_obj, greedy_result$subset, maximize = TRUE)
cat("Exchange refinement completed\n")
cat("Improved subset size:", length(exchange_result$subset), "\n")
cat("Iterations:", exchange_result$iterations, "\n")
cat("Improvements found:", length(exchange_result$improvements), "\n")
cat("Converged:", exchange_result$converged, "\n")

# Test 6: Complete algorithm
cat("\n7. Testing complete algorithm...\n")
source("R/subset_greedy.R")
source("R/subset_exchange.R")

complete_result <- run_complete_algorithm(dist_obj, subset_size = 5, maximize = TRUE)
cat("Complete algorithm completed\n")
cat("Final subset:", complete_result$final_subset_names, "\n")
cat("Final metrics:\n")
print(complete_result$final_metrics)

# Test 7: Random subset sampling
cat("\n8. Testing random subset sampling...\n")
source("R/subset_random.R")

random_subsets <- sample_random_subsets(dist_obj, subset_size = 5, n_reps = 10)
cat("Generated", length(random_subsets), "random subsets\n")

# Test 8: Result 1 analysis
cat("\n9. Testing Result 1 analysis...\n")
source("R/result1_analysis.R")

# Create distance objects for all trees
dist_objs <- list(
  small_balanced = dist_obj
)

result1 <- run_result1_analysis(dist_objs, test_cfg)
cat("Result 1 analysis completed\n")
cat("Number of trees analyzed:", length(result1$results), "\n")

# Test 9: Result 2 analysis
cat("\n10. Testing Result 2 analysis...\n")
source("R/result2_analysis.R")

result2 <- run_result2_analysis(dist_objs, test_cfg)
cat("Result 2 analysis completed\n")

# Test 10: Plotting functions
cat("\n11. Testing plotting functions...\n")
source("R/plotting.R")

# Test basic plotting
library(ggplot2)
test_data <- data.frame(
  Algorithm = c("A", "B", "C"),
  MinPD = c(10, 9, 11),
  MeanPD = c(20, 22, 19),
  MeanNND = c(5, 6, 4)
)

p <- plot_algorithm_comparison(test_data, title = "Test Plot")
cat("Plot created successfully\n")

# Test 11: Utility functions
cat("\n12. Testing utility functions...\n")
source("R/utils_io.R")

# Create test directories
dirs <- create_output_dirs("test_output")
cat("Created directories:", names(dirs), "\n")

# Test saving and loading
test_df <- data.frame(x = 1:5, y = letters[1:5])
saved_path <- save_csv(test_df, "test_data", "test_output/tables")
loaded_df <- load_csv("test_data", "test_output/tables")
cat("Data frame saved and loaded successfully\n")
cat("Original and loaded identical:", identical(test_df, loaded_df), "\n")

# Clean up
unlink("test_output", recursive = TRUE)
cat("Cleaned up test directories\n")

# Summary
cat(paste0("\n", strrep("=", 80), "\n"))
cat("IMPLEMENTATION TEST COMPLETED SUCCESSFULLY\n")
cat(paste0(strrep("=", 80), "\n"))

cat("\nAll modules tested:\n")
cat("✓ Configuration\n")
cat("✓ Tree generation\n")
cat("✓ Distance metrics\n")
cat("✓ Objective comparison\n")
cat("✓ Greedy algorithm\n")
cat("✓ Exchange refinement\n")
cat("✓ Complete algorithm\n")
cat("✓ Random subset sampling\n")
cat("✓ Result 1 analysis\n")
cat("✓ Result 2 analysis\n")
cat("✓ Plotting functions\n")
cat("✓ Utility functions\n")

cat("\nThe implementation appears to be working correctly.\n")
cat("For a full test, run: Rscript main.R test\n")
