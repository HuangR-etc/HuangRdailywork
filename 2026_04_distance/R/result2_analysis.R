# Result 2 analysis: Test design components
# This module implements the analysis for Result 2

#' Run multiple algorithms for comparison
#'
#' @param dist_obj Distance object
#' @param subset_size Subset size
#' @param n_greedy_starts Number of random starts for greedy phase
#' @return Results from all algorithms
run_all_algorithms <- function(dist_obj, subset_size, n_greedy_starts = 1) {
  cat("Running all algorithms for comparison...\n")
  
  results <- list()
  
  # Algorithm A: Main method (MinPD > MeanPD > MeanNND, greedy + exchange)
  cat("  Algorithm A: Main method (lexicographic, greedy + exchange)...\n")
  results$A_main <- run_complete_algorithm(dist_obj, subset_size, maximize = TRUE,
                                           n_greedy_starts = n_greedy_starts)
  
  # Algorithm B: Greedy only (MinPD > MeanPD > MeanNND, greedy only)
  cat("  Algorithm B: Greedy only (lexicographic, greedy only)...\n")
  greedy_result <- build_subset_greedy(dist_obj, subset_size, maximize = TRUE)
  results$B_greedy_only <- list(
    final_subset = greedy_result$subset,
    final_subset_names = greedy_result$subset_names,
    final_metrics = greedy_result$metrics,
    algorithm = greedy_result$algorithm
  )
  
  # Algorithm C: MeanPD only (maximize MeanPD, greedy + exchange)
  cat("  Algorithm C: MeanPD only (greedy + exchange)...\n")
  results$C_meanpd_only <- run_complete_algorithm(dist_obj, subset_size, maximize = TRUE,
                                                  single_objective = "MeanPD",
                                                  n_greedy_starts = n_greedy_starts)
  
  # Algorithm D: MinPD only (maximize MinPD, greedy + exchange)
  cat("  Algorithm D: MinPD only (greedy + exchange)...\n")
  results$D_minpd_only <- run_complete_algorithm(dist_obj, subset_size, maximize = TRUE,
                                                 single_objective = "MinPD",
                                                 n_greedy_starts = n_greedy_starts)
  
  # Algorithm E: MeanNND only (maximize MeanNND, greedy + exchange)
  cat("  Algorithm E: MeanNND only (greedy + exchange)...\n")
  results$E_meannnd_only <- run_complete_algorithm(dist_obj, subset_size, maximize = TRUE,
                                                   single_objective = "MeanNND",
                                                   n_greedy_starts = n_greedy_starts)
  
  return(results)
}

#' Compare algorithm results
#'
#' @param algorithm_results Results from run_all_algorithms
#' @return Comparison metrics
compare_algorithms <- function(algorithm_results) {
  cat("Comparing algorithm results...\n")
  
  # Extract metrics from each algorithm
  metrics_df <- data.frame()
  
  for (algo_name in names(algorithm_results)) {
    result <- algorithm_results[[algo_name]]
    
    metrics_df <- rbind(metrics_df, data.frame(
      Algorithm = algo_name,
      MinPD = result$final_metrics$MinPD,
      MeanPD = result$final_metrics$MeanPD,
      MeanNND = result$final_metrics$MeanNND,
      Subset_Size = length(result$final_subset),
      stringsAsFactors = FALSE
    ))
  }
  
  # Calculate pairwise comparisons
  comparisons <- list()
  
  # A vs B: Effect of exchange refinement
  a_metrics <- algorithm_results$A_main$final_metrics
  b_metrics <- algorithm_results$B_greedy_only$final_metrics
  
  comparisons$A_vs_B <- list(
    MinPD_diff = a_metrics$MinPD - b_metrics$MinPD,
    MeanPD_diff = a_metrics$MeanPD - b_metrics$MeanPD,
    MeanNND_diff = a_metrics$MeanNND - b_metrics$MeanNND,
    A_better_MinPD = is_better_lexico_max(a_metrics, b_metrics)
  )
  
  # A vs C/D/E: Effect of multi-criterion vs single objective
  for (algo in c("C_meanpd_only", "D_minpd_only", "E_meannnd_only")) {
    other_metrics <- algorithm_results[[algo]]$final_metrics
    comp_name <- paste0("A_vs_", algo)
    
    comparisons[[comp_name]] <- list(
      MinPD_diff = a_metrics$MinPD - other_metrics$MinPD,
      MeanPD_diff = a_metrics$MeanPD - other_metrics$MeanPD,
      MeanNND_diff = a_metrics$MeanNND - other_metrics$MeanNND,
      A_better = is_better_lexico_max(a_metrics, other_metrics)
    )
  }
  
  # Compare single objective algorithms with each other
  comparisons$C_vs_D <- list(
    MinPD_diff = algorithm_results$C_meanpd_only$final_metrics$MinPD - 
                 algorithm_results$D_minpd_only$final_metrics$MinPD,
    MeanPD_diff = algorithm_results$C_meanpd_only$final_metrics$MeanPD - 
                  algorithm_results$D_minpd_only$final_metrics$MeanPD,
    MeanNND_diff = algorithm_results$C_meanpd_only$final_metrics$MeanNND - 
                   algorithm_results$D_minpd_only$final_metrics$MeanNND
  )
  
  comparisons$C_vs_E <- list(
    MinPD_diff = algorithm_results$C_meanpd_only$final_metrics$MinPD - 
                 algorithm_results$E_meannnd_only$final_metrics$MinPD,
    MeanPD_diff = algorithm_results$C_meanpd_only$final_metrics$MeanPD - 
                  algorithm_results$E_meannnd_only$final_metrics$MeanPD,
    MeanNND_diff = algorithm_results$C_meanpd_only$final_metrics$MeanNND - 
                   algorithm_results$E_meannnd_only$final_metrics$MeanNND
  )
  
  comparisons$D_vs_E <- list(
    MinPD_diff = algorithm_results$D_minpd_only$final_metrics$MinPD - 
                 algorithm_results$E_meannnd_only$final_metrics$MinPD,
    MeanPD_diff = algorithm_results$D_minpd_only$final_metrics$MeanPD - 
                  algorithm_results$E_meannnd_only$final_metrics$MeanPD,
    MeanNND_diff = algorithm_results$D_minpd_only$final_metrics$MeanNND - 
                   algorithm_results$E_meannnd_only$final_metrics$MeanNND
  )
  
  return(list(
    metrics = metrics_df,
    comparisons = comparisons
  ))
}

#' Run Result 2 analysis for a single tree
#'
#' @param dist_obj Distance object
#' @param subset_size Subset size
#' @param n_greedy_starts Number of random starts for greedy phase
#' @return Comprehensive Result 2 analysis results
run_result2_for_tree <- function(dist_obj, subset_size, n_greedy_starts = 1) {
  
  tree_name <- ifelse("tree_name" %in% names(dist_obj), 
                     dist_obj$tree_name, 
                     paste0("tree_", length(dist_obj$tip_labels), "tips"))
  
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat("Result 2 analysis for", tree_name, "\n")
  cat("Subset size:", subset_size, "\n")
  cat(paste0(strrep("=", 60), "\n\n"))
  
  # Step 1: Run all algorithms
  cat("Step 1: Running all algorithms...\n")
  algorithm_results <- run_all_algorithms(dist_obj, subset_size, n_greedy_starts)
  
  # Step 2: Compare algorithms
  cat("\nStep 2: Comparing algorithms...\n")
  comparison <- compare_algorithms(algorithm_results)
  
  # Step 3: Create summary
  summary <- list(
    tree_name = tree_name,
    tree_n_tips = length(dist_obj$tip_labels),
    subset_size = subset_size,
    algorithm_metrics = comparison$metrics,
    comparisons = comparison$comparisons
  )
  
  # Step 4: Print results
  cat("\nStep 3: Algorithm performance\n")
  cat(paste0(strrep("-", 60), "\n"))
  cat("Algorithm        | MinPD   | MeanPD  | MeanNND | Better?\n")
  cat(paste0(strrep("-", 60), "\n"))
  
  for (i in 1:nrow(comparison$metrics)) {
    row <- comparison$metrics[i, ]
    algo_name <- row$Algorithm
    
    # Determine if this algorithm is the main algorithm (A)
    is_main <- (algo_name == "A_main")
    
    cat(sprintf("%-16s | %7.3f | %7.3f | %7.3f | %s\n",
                algo_name, row$MinPD, row$MeanPD, row$MeanNND,
                ifelse(is_main, "MAIN", "")))
  }
  
  cat(paste0(strrep("-", 60), "\n"))
  
  # Print key comparisons
  cat("\nStep 4: Key comparisons\n")
  cat(paste0(strrep("-", 60), "\n"))
  
  # A vs B: Exchange refinement
  a_vs_b <- comparison$comparisons$A_vs_B
  cat("A (main) vs B (greedy only):\n")
  cat("  MinPD difference:  ", sprintf("%+7.3f", a_vs_b$MinPD_diff), 
      ifelse(a_vs_b$MinPD_diff > 0, "(improvement)", "(worse)"), "\n")
  cat("  MeanPD difference: ", sprintf("%+7.3f", a_vs_b$MeanPD_diff), 
      ifelse(a_vs_b$MeanPD_diff > 0, "(improvement)", "(worse)"), "\n")
  cat("  MeanNND difference:", sprintf("%+7.3f", a_vs_b$MeanNND_diff), 
      ifelse(a_vs_b$MeanNND_diff > 0, "(improvement)", "(worse)"), "\n")
  cat("  A is lexicographically better than B:", a_vs_b$A_better_MinPD, "\n")
  
  # A vs single objective algorithms
  cat("\nA (main) vs single objective algorithms:\n")
  for (algo in c("C_meanpd_only", "D_minpd_only", "E_meannnd_only")) {
    comp <- comparison$comparisons[[paste0("A_vs_", algo)]]
    cat("  vs", algo, ":", comp$A_better, "\n")
  }
  
  cat(paste0(strrep("-", 60), "\n"))
  
  # Conclusions
  cat("\nStep 5: Conclusions\n")
  cat(paste0(strrep("-", 60), "\n"))
  
  # Check if exchange refinement helps
  if (a_vs_b$A_better_MinPD) {
    cat("1. Exchange refinement IMPROVES solution quality (A better than B)\n")
  } else {
    cat("1. Exchange refinement does NOT improve solution quality\n")
  }
  
  # Check if multi-criterion is better than single objectives
  better_count <- 0
  for (algo in c("C_meanpd_only", "D_minpd_only", "E_meannnd_only")) {
    comp <- comparison$comparisons[[paste0("A_vs_", algo)]]
    if (comp$A_better) better_count <- better_count + 1
  }
  
  if (better_count == 3) {
    cat("2. Multi-criterion optimization is BETTER than all single objectives\n")
  } else if (better_count > 0) {
    cat("2. Multi-criterion optimization is better than", better_count, "out of 3 single objectives\n")
  } else {
    cat("2. Multi-criterion optimization is NOT better than single objectives\n")
  }
  
  # Check trade-offs between single objectives
  cat("3. Single objective trade-offs:\n")
  cat("   - Maximizing MeanPD alone gives highest MeanPD:", 
      max(comparison$metrics$MeanPD[comparison$metrics$Algorithm != "A_main"]), "\n")
  cat("   - Maximizing MinPD alone gives highest MinPD:", 
      max(comparison$metrics$MinPD[comparison$metrics$Algorithm != "A_main"]), "\n")
  cat("   - Maximizing MeanNND alone gives highest MeanNND:", 
      max(comparison$metrics$MeanNND[comparison$metrics$Algorithm != "A_main"]), "\n")
  
  return(list(
    algorithm_results = algorithm_results,
    comparison = comparison,
    summary = summary
  ))
}

#' Run Result 2 analysis for all trees
#'
#' @param dist_objs List of distance objects for all trees
#' @param cfg Configuration list
#' @return Comprehensive Result 2 results for all trees
run_result2_analysis <- function(dist_objs, cfg) {
  cat(paste0("\n", strrep("=", 70), "\n"))
  cat("RESULT 2: Testing design components\n")
  cat(paste0(strrep("=", 70), "\n"))
  
  results <- list()
  
  # Process large balanced tree
  if ("large_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing large balanced tree (n =", cfg$large_n, ")\n")
    result_lb <- run_result2_for_tree(
      dist_objs$large_balanced,
      subset_size = cfg$subset_large,
      n_greedy_starts = 1
    )
    results$large_balanced <- result_lb
  }
  
  # Process large ladder tree
  if ("large_ladder" %in% names(dist_objs)) {
    cat("\n>>> Processing large ladder tree (n =", cfg$large_n, ")\n")
    result_ll <- run_result2_for_tree(
      dist_objs$large_ladder,
      subset_size = cfg$subset_large,
      n_greedy_starts = 1
    )
    results$large_ladder <- result_ll
  }
  
  # Process small balanced tree
  if ("small_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing small balanced tree (n =", cfg$small_n, ")\n")
    result_sb <- run_result2_for_tree(
      dist_objs$small_balanced,
      subset_size = cfg$subset_small,
      n_greedy_starts = 1
    )
    results$small_balanced <- result_sb
  }
  
  # Create overall summary
  overall_summary <- data.frame()
  
  for (tree_name in names(results)) {
    result <- results[[tree_name]]
    metrics <- result$comparison$metrics
    
    for (i in 1:nrow(metrics)) {
      overall_summary <- rbind(overall_summary, data.frame(
        Tree = tree_name,
        Algorithm = metrics$Algorithm[i],
        MinPD = metrics$MinPD[i],
        MeanPD = metrics$MeanPD[i],
        MeanNND = metrics$MeanNND[i],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat(paste0("\n", strrep("=", 70), "\n"))
  cat("OVERALL SUMMARY\n")
  cat(paste0(strrep("=", 70), "\n"))
  print(overall_summary)
  
  return(list(
    results = results,
    overall_summary = overall_summary
  ))
}

#' Save Result 2 results to files
#'
#' @param result2_results Results from run_result2_analysis
#' @param cfg Configuration list
save_result2 <- function(result2_results, cfg) {
  cat("\nSaving Result 2 results...\n")
  
  # Create output directory if it doesn't exist
  output_dir <- cfg$tables_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save overall summary
  summary_file <- file.path(output_dir, "result2_summary.csv")
  write.csv(result2_results$overall_summary, summary_file, row.names = FALSE)
  cat("  Saved overall summary to:", summary_file, "\n")
  
  # Save detailed results for each tree
  for (tree_name in names(result2_results$results)) {
    result <- result2_results$results[[tree_name]]
    
    # Save algorithm metrics
    metrics_file <- file.path(output_dir, paste0("result2_metrics_", tree_name, ".csv"))
    write.csv(result$comparison$metrics, metrics_file, row.names = FALSE)
    
    # Save comparison results
    comp_file <- file.path(output_dir, paste0("result2_comparisons_", tree_name, ".csv"))
    
    # Convert comparisons to data frame
    comp_list <- result$comparison$comparisons
    comp_df <- data.frame()
    
    for (comp_name in names(comp_list)) {
      comp <- comp_list[[comp_name]]
      comp_df <- rbind(comp_df, data.frame(
        Comparison = comp_name,
        MinPD_diff = ifelse("MinPD_diff" %in% names(comp), comp$MinPD_diff, NA),
        MeanPD_diff = ifelse("MeanPD_diff" %in% names(comp), comp$MeanPD_diff, NA),
        MeanNND_diff = ifelse("MeanNND_diff" %in% names(comp), comp$MeanNND_diff, NA),
        A_better = ifelse("A_better" %in% names(comp), comp$A_better, 
                         ifelse("A_better_MinPD" %in% names(comp), comp$A_better_MinPD, NA)),
        stringsAsFactors = FALSE
      ))
    }
    
    write.csv(comp_df, comp_file, row.names = FALSE)
  }
  
  cat("Result 2 results saved to:", output_dir, "\n")
}

#' Test Result 2 analysis
test_result2 <- function() {
  # Load required libraries and functions
  library(ape)
  source("config/analysis_config.R")
  source("distance_metrics.R")
  source("objective_compare.R")
  source("subset_greedy.R")
  source("subset_exchange.R")
  
  # Create a test configuration
  test_cfg <- list(
    large_n = 30,  # Smaller for testing
    small_n = 15,
    subset_large = 5,
    subset_small = 3
  )
  
  # Create a test tree
  cat("Creating test tree...\n")
  test_tree <- rtree(test_cfg$large_n)
  test_tree$tip.label <- paste0("sp", 1:test_cfg$large_n)
  dist_obj <- create_distance_object(test_tree)
  dist_obj$tree_name <- "test_tree"
  
  # Run Result 2 analysis
  cat("\nRunning Result 2 analysis...\n")
  result2 <- run_result2_for_tree(dist_obj, subset_size = test_cfg$subset_large)
  
  # Save results
  test_cfg$tables_dir <- "test_output"
  if (!dir.exists(test_cfg$tables_dir)) {
    dir.create(test_cfg$tables_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create a simple result structure for saving
  simple_result <- list(
    results = list(test_tree = result2),
    overall_summary = result2$comparison$metrics
  )
  
  save_result2(simple_result, test_cfg)
  
  return(result2)
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_result <- test_result2()
  cat("\nTest completed successfully.\n")
}
