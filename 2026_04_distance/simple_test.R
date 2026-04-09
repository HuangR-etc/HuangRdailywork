# Simple test to verify core functionality
cat("Simple test of core functionality\n")

# Load configuration
source("config/analysis_config.R")

# Create test configuration
test_cfg <- cfg
test_cfg$large_n <- 32
test_cfg$small_n <- 16
test_cfg$subset_large <- 5
test_cfg$subset_small <- 3
test_cfg$null_reps_large <- 10
test_cfg$null_reps_small <- 10
test_cfg$trait_reps <- 2
test_cfg$random_subset_reps_for_trait <- 5
test_cfg$exhaustive_small <- FALSE

# Set seed
set.seed(test_cfg$seed)

# Load core modules
cat("Loading core modules...\n")
source("R/tree_generators.R")
source("R/distance_metrics.R")
source("R/objective_compare.R")
source("R/subset_greedy.R")
source("R/subset_exchange.R")
source("R/subset_random.R")

# Test 1: Tree generation
cat("\n1. Testing tree generation...\n")
trees <- generate_all_trees(test_cfg)
cat("Generated trees:", names(trees), "\n")

# Test 2: Distance metrics
cat("\n2. Testing distance metrics...\n")
dist_obj <- create_distance_object(trees$small_balanced)
cat("Distance matrix dimensions:", dim(dist_obj$dist_mat), "\n")

# Test 3: Greedy algorithm
cat("\n3. Testing greedy algorithm...\n")
greedy_result <- build_subset_greedy(dist_obj, subset_size = 5, maximize = TRUE)
cat("Greedy subset size:", length(greedy_result$subset), "\n")
cat("Greedy metrics: MinPD =", greedy_result$metrics$MinPD, 
    "MeanPD =", greedy_result$metrics$MeanPD, 
    "MeanNND =", greedy_result$metrics$MeanNND, "\n")

# Test 4: Exchange refinement
cat("\n4. Testing exchange refinement...\n")
exchange_result <- refine_subset_exchange(dist_obj, greedy_result$subset, maximize = TRUE)
cat("Exchange subset size:", length(exchange_result$subset), "\n")
cat("Exchange metrics: MinPD =", exchange_result$metrics$MinPD, 
    "MeanPD =", exchange_result$metrics$MeanPD, 
    "MeanNND =", exchange_result$metrics$MeanNND, "\n")

# Test 5: Random subsets
cat("\n5. Testing random subsets...\n")
random_subsets <- sample_random_subsets(dist_obj, subset_size = 5, n_reps = 5)
cat("Generated", length(random_subsets), "random subsets\n")

# Test 6: Objective comparison
cat("\n6. Testing objective comparison...\n")
metrics1 <- list(MinPD = 10, MeanPD = 20, MeanNND = 5)
metrics2 <- list(MinPD = 9, MeanPD = 25, MeanNND = 6)
cat("Is metrics1 better than metrics2?", is_better_lexico_max(metrics1, metrics2), "\n")
cat("Is metrics2 better than metrics1?", is_better_lexico_max(metrics2, metrics1), "\n")

# Test 7: Complete algorithm
cat("\n7. Testing complete algorithm...\n")
complete_result <- run_complete_algorithm(dist_obj, subset_size = 5, maximize = TRUE)
cat("Complete algorithm completed\n")
cat("Final subset size:", length(complete_result$final_subset), "\n")
cat("Final metrics: MinPD =", complete_result$final_metrics$MinPD, 
    "MeanPD =", complete_result$final_metrics$MeanPD, 
    "MeanNND =", complete_result$final_metrics$MeanNND, "\n")

cat("\nAll core functionality tests passed!\n")
