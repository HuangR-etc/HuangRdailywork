# Test script to verify signal metrics calculation fix
library(ape)
library(picante)
library(phytools)

# Source the fixed functions
source("R/signal_metrics.R")
source("R/trait_simulation.R")

cat("Testing fixed signal metrics calculation...\n")

# Create a test tree
set.seed(123)
test_tree <- rtree(20)
test_tree$tip.label <- paste0("sp", 1:20)

# Create test subsets
subsets <- list(
  c("sp1", "sp3", "sp5", "sp7", "sp9"),
  c("sp2", "sp4", "sp6", "sp8", "sp10"),
  c("sp1", "sp5", "sp10", "sp15", "sp20")
)
subset_names <- c("subset1", "subset2", "subset3")

# Test 1: Test calculate_subsets_signal function
cat("\n1. Testing calculate_subsets_signal function:\n")
# Simulate trait data under BM
trait_values <- fastBM(test_tree, sig2 = 1.0)
names(trait_values) <- test_tree$tip.label

subset_signal <- calculate_subsets_signal(test_tree, subsets, subset_names, 
                                         trait_values, metrics = c("K", "lambda"))
print(subset_signal)

# Check if any values are NA
if (any(is.na(subset_signal$K)) || any(is.na(subset_signal$Lambda))) {
  cat("WARNING: Some signal values are NA!\n")
} else {
  cat("SUCCESS: All signal values calculated successfully!\n")
}

# Test 2: Test simulate_traits_for_subsets and analyze_simulated_traits_signal
cat("\n2. Testing simulate_traits_for_subsets and analyze_simulated_traits_signal:\n")
sim_results <- simulate_traits_for_subsets(test_tree, subsets, subset_names, 
                                          n_reps = 5, ou_alpha_values = c(0.2, 1))

# Check if subsets are stored in metadata
if (!is.null(sim_results$metadata$subsets)) {
  cat("SUCCESS: Subsets are stored in metadata\n")
} else {
  cat("ERROR: Subsets are NOT stored in metadata\n")
}

# Analyze signal
signal_results <- analyze_simulated_traits_signal(sim_results, test_tree, n_reps_to_analyze = 3)

cat("Signal results (first few rows):\n")
print(head(signal_results))

# Check for NA values
na_count <- sum(is.na(signal_results$K) | is.na(signal_results$Lambda))
total_count <- nrow(signal_results)

if (na_count > 0) {
  cat(paste("WARNING:", na_count, "out of", total_count, "values are NA\n"))
  
  # Check which rows have NA
  na_rows <- signal_results[is.na(signal_results$K) | is.na(signal_results$Lambda), ]
  cat("Rows with NA values:\n")
  print(na_rows)
} else {
  cat(paste("SUCCESS: All", total_count, "signal values calculated successfully!\n"))
}

# Test 3: Test individual functions
cat("\n3. Testing individual calculation functions:\n")

# Test calc_blombergs_K
test_subset <- subsets[[1]]
test_subset_tree <- keep.tip(test_tree, test_subset)
test_trait_values <- trait_values[test_subset]

K_value <- calc_blombergs_K(test_trait_values, test_subset_tree)
cat(paste("  Blomberg's K for subset1:", K_value, "\n"))

# Test calc_pagels_lambda
lambda_value <- calc_pagels_lambda(test_trait_values, test_subset_tree)
cat(paste("  Pagel's lambda for subset1:", lambda_value, "\n"))

if (!is.na(K_value) && !is.na(lambda_value)) {
  cat("  SUCCESS: Individual functions work correctly\n")
} else {
  cat("  WARNING: Individual functions returned NA\n")
}

cat("\nTest completed.\n")
