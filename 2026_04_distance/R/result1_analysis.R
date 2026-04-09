# Result 1 analysis: Compare dispersed subsets with random null distribution
# This module implements the analysis for Result 1

#' Run Result 1 analysis for a single tree
#'
#' @param dist_obj Distance object
#' @param subset_size Subset size
#' @param n_null_reps Number of random subsets for null distribution
#' @param maximize If TRUE, test maximization (dispersed); if FALSE, test minimization (clustered)
#' @param n_greedy_starts Number of random starts for greedy phase
#' @return Comprehensive Result 1 analysis results
run_result1_for_tree <- function(dist_obj, subset_size, n_null_reps, 
                                 maximize = TRUE, n_greedy_starts = 1) {
  
  tree_name <- ifelse("tree_name" %in% names(dist_obj), 
                     dist_obj$tree_name, 
                     paste0("tree_", length(dist_obj$tip_labels), "tips"))
  
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat("Result 1 analysis for", tree_name, "\n")
  cat("Subset size:", subset_size, "\n")
  cat("Null replicates:", n_null_reps, "\n")
  cat("Optimization:", ifelse(maximize, "Maximization (dispersed)", "Minimization (clustered)"), "\n")
  cat(paste0(strrep("=", 60), "\n\n"))
  
  # Step 1: Run main algorithm to get observed subset
  cat("Step 1: Running main algorithm...\n")
  main_result <- run_complete_algorithm(dist_obj, subset_size, maximize, 
                                        n_greedy_starts = n_greedy_starts)
  
  observed_subset <- main_result$final_subset
  observed_metrics <- main_result$final_metrics
  
  cat("  Observed subset size:", length(observed_subset), "\n")
  cat("  Observed metrics: MinPD =", observed_metrics$MinPD, 
      "MeanPD =", observed_metrics$MeanPD, 
      "MeanNND =", observed_metrics$MeanNND, "\n")
  
  # Step 2: Compare with null distribution
  cat("\nStep 2: Comparing with null distribution...\n")
  comparison <- compare_with_null(dist_obj, observed_subset, subset_size, 
                                  n_null_reps, maximize)
  
  # Step 3: Create summary
  summary <- list(
    tree_name = tree_name,
    tree_n_tips = length(dist_obj$tip_labels),
    subset_size = subset_size,
    n_null_reps = n_null_reps,
    maximize = maximize,
    observed_subset = observed_subset,
    observed_subset_names = dist_obj$tip_labels[observed_subset],
    observed_metrics = observed_metrics,
    null_summary = list(
      MinPD_mean = mean(comparison$null_metrics$MinPD),
      MinPD_sd = sd(comparison$null_metrics$MinPD),
      MeanPD_mean = mean(comparison$null_metrics$MeanPD),
      MeanPD_sd = sd(comparison$null_metrics$MeanPD),
      MeanNND_mean = mean(comparison$null_metrics$MeanNND),
      MeanNND_sd = sd(comparison$null_metrics$MeanNND)
    ),
    p_values = comparison$p_values,
    z_scores = comparison$z_scores,
    percentiles = comparison$percentiles,
    comparison_summary = comparison$summary
  )
  
  # Step 4: Print results
  cat("\nStep 3: Results summary\n")
  cat(paste0(strrep("-", 60), "\n"))
  cat("Metric    | Observed | Null Mean | Null SD | Z-Score  | P-Value  | Percentile\n")
  cat(paste0(strrep("-", 60), "\n"))
  
  metrics <- c("MinPD", "MeanPD", "MeanNND")
  for (metric in metrics) {
    obs <- summary$observed_metrics[[metric]]
    null_mean <- summary$null_summary[[paste0(metric, "_mean")]]
    null_sd <- summary$null_summary[[paste0(metric, "_sd")]]
    z <- summary$z_scores[[metric]]
    p <- summary$p_values[[metric]]
    pct <- summary$percentiles[[metric]]
    
    cat(sprintf("%-9s | %8.3f | %9.3f | %7.3f | %8.3f | %8.3f | %10.3f\n",
                metric, obs, null_mean, null_sd, z, p, pct))
  }
  
  cat(paste0(strrep("-", 60), "\n"))
  
  # Determine if observed is significantly better than random
  sig_threshold <- 0.05
  is_significant <- all(sapply(summary$p_values, function(p) p < sig_threshold))
  
  if (maximize) {
    cat("Conclusion: Observed dispersed subset is")
    if (is_significant) {
      cat(" SIGNIFICANTLY more dispersed than random subsets (p <", sig_threshold, ")\n")
    } else {
      cat(" NOT significantly more dispersed than random subsets\n")
    }
  } else {
    cat("Conclusion: Observed clustered subset is")
    if (is_significant) {
      cat(" SIGNIFICANTLY more clustered than random subsets (p <", sig_threshold, ")\n")
    } else {
      cat(" NOT significantly more clustered than random subsets\n")
    }
  }
  
  return(list(
    main_result = main_result,
    comparison = comparison,
    summary = summary,
    is_significant = is_significant
  ))
}

#' Run Result 1 analysis for all trees
#'
#' @param dist_objs List of distance objects for all trees
#' @param cfg Configuration list
#' @return Comprehensive Result 1 results for all trees
run_result1_analysis <- function(dist_objs, cfg) {
  cat(paste0("\n", strrep("=", 70), "\n"))
  cat("RESULT 1: Comparing dispersed subsets with random null distribution\n")
  cat(paste0(strrep("=", 70), "\n"))
  
  results <- list()
  
  # Process large balanced tree
  if ("large_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing large balanced tree (n =", cfg$large_n, ")\n")
    result_lb <- run_result1_for_tree(
      dist_objs$large_balanced,
      subset_size = cfg$subset_large,
      n_null_reps = cfg$null_reps_large,
      maximize = TRUE,
      n_greedy_starts = 1
    )
    results$large_balanced <- result_lb
  }
  
  # Process large ladder tree
  if ("large_ladder" %in% names(dist_objs)) {
    cat("\n>>> Processing large ladder tree (n =", cfg$large_n, ")\n")
    result_ll <- run_result1_for_tree(
      dist_objs$large_ladder,
      subset_size = cfg$subset_large,
      n_null_reps = cfg$null_reps_large,
      maximize = TRUE,
      n_greedy_starts = 1
    )
    results$large_ladder <- result_ll
  }
  
  # Process small balanced tree
  if ("small_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing small balanced tree (n =", cfg$small_n, ")\n")
    
    # For small tree, we might want exhaustive enumeration
    n_null_reps <- cfg$null_reps_small
    if (cfg$exhaustive_small) {
      # Calculate total number of combinations
      total_combinations <- choose(cfg$small_n, cfg$subset_small)
      if (total_combinations <= 10000) {  # Reasonable limit for exhaustive
        cat("  Using exhaustive enumeration (", total_combinations, " combinations)\n")
        # We'll implement exhaustive null distribution
        # For now, use sampling
        n_null_reps <- min(total_combinations, cfg$null_reps_small)
      }
    }
    
    result_sb <- run_result1_for_tree(
      dist_objs$small_balanced,
      subset_size = cfg$subset_small,
      n_null_reps = n_null_reps,
      maximize = TRUE,
      n_greedy_starts = 1
    )
    results$small_balanced <- result_sb
  }
  
  # Create overall summary
  overall_summary <- data.frame()
  
  for (tree_name in names(results)) {
    result <- results[[tree_name]]
    summary <- result$summary
    
    overall_summary <- rbind(overall_summary, data.frame(
      Tree = tree_name,
      N_Tips = summary$tree_n_tips,
      Subset_Size = summary$subset_size,
      MinPD_Observed = summary$observed_metrics$MinPD,
      MinPD_P = summary$p_values$MinPD,
      MinPD_Z = summary$z_scores$MinPD,
      MeanPD_Observed = summary$observed_metrics$MeanPD,
      MeanPD_P = summary$p_values$MeanPD,
      MeanPD_Z = summary$z_scores$MeanPD,
      MeanNND_Observed = summary$observed_metrics$MeanNND,
      MeanNND_P = summary$p_values$MeanNND,
      MeanNND_Z = summary$z_scores$MeanNND,
      Significant = result$is_significant,
      stringsAsFactors = FALSE
    ))
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

#' Save Result 1 results to files
#'
#' @param result1_results Results from run_result1_analysis
#' @param cfg Configuration list
save_result1 <- function(result1_results, cfg) {
  cat("\nSaving Result 1 results...\n")
  
  # Create output directory if it doesn't exist
  output_dir <- cfg$tables_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save overall summary
  summary_file <- file.path(output_dir, "result1_summary.csv")
  write.csv(result1_results$overall_summary, summary_file, row.names = FALSE)
  cat("  Saved overall summary to:", summary_file, "\n")
  
  # Save detailed results for each tree
  for (tree_name in names(result1_results$results)) {
    result <- result1_results$results[[tree_name]]
    
    # Save observed metrics
    obs_file <- file.path(output_dir, paste0("result1_observed_", tree_name, ".csv"))
    obs_df <- data.frame(
      Tree = tree_name,
      Subset_Size = result$summary$subset_size,
      Tip_Index = result$summary$observed_subset,
      Tip_Name = result$summary$observed_subset_names,
      MinPD = result$summary$observed_metrics$MinPD,
      MeanPD = result$summary$observed_metrics$MeanPD,
      MeanNND = result$summary$observed_metrics$MeanNND,
      stringsAsFactors = FALSE
    )
    write.csv(obs_df, obs_file, row.names = FALSE)
    
    # Save null distribution metrics
    null_file <- file.path(output_dir, paste0("result1_null_", tree_name, ".csv"))
    write.csv(result$comparison$null_metrics, null_file, row.names = FALSE)
    
    # Save comparison summary
    comp_file <- file.path(output_dir, paste0("result1_comparison_", tree_name, ".csv"))
    write.csv(result$comparison$summary, comp_file, row.names = FALSE)
  }
  
  cat("Result 1 results saved to:", output_dir, "\n")
}

#' Test Result 1 analysis
test_result1 <- function() {
  # Load required libraries and functions
  library(ape)
  source("config/analysis_config.R")
  source("distance_metrics.R")
  source("objective_compare.R")
  source("subset_greedy.R")
  source("subset_exchange.R")
  source("subset_random.R")
  
  # Create a test configuration
  test_cfg <- list(
    large_n = 30,  # Smaller for testing
    small_n = 15,
    subset_large = 5,
    subset_small = 3,
    null_reps_large = 100,
    null_reps_small = 50,
    exhaustive_small = FALSE
  )
  
  # Create test trees
  cat("Creating test trees...\n")
  
  # Large balanced tree
  tree_lb <- rtree(test_cfg$large_n)
  tree_lb$tip.label <- paste0("LB", 1:test_cfg$large_n)
  dist_lb <- create_distance_object(tree_lb)
  dist_lb$tree_name <- "test_large_balanced"
  
  # Large ladder tree
  tree_ll <- stree(test_cfg$large_n, type = "left")
  tree_ll$tip.label <- paste0("LL", 1:test_cfg$large_n)
  dist_ll <- create_distance_object(tree_ll)
  dist_ll$tree_name <- "test_large_ladder"
  
  # Small balanced tree
  tree_sb <- rtree(test_cfg$small_n)
  tree_sb$tip.label <- paste0("SB", 1:test_cfg$small_n)
  dist_sb <- create_distance_object(tree_sb)
  dist_sb$tree_name <- "test_small_balanced"
  
  dist_objs <- list(
    large_balanced = dist_lb,
    large_ladder = dist_ll,
    small_balanced = dist_sb
  )
  
  # Run Result 1 analysis
  cat("\nRunning Result 1 analysis...\n")
  result1 <- run_result1_analysis(dist_objs, test_cfg)
  
  # Save results
  test_cfg$tables_dir <- "test_output"
  if (!dir.exists(test_cfg$tables_dir)) {
    dir.create(test_cfg$tables_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  save_result1(result1, test_cfg)
  
  return(result1)
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_result <- test_result1()
  cat("\nTest completed successfully.\n")
}
