# Result 4 analysis: Compare heuristic with exact optimum
# This module implements the analysis for Result 4

#' Run Result 4 analysis for small tree
#'
#' @param dist_obj Distance object (for small tree)
#' @param subset_size Subset size
#' @param maximize If TRUE, compare for maximization (default: TRUE)
#' @param n_greedy_starts Number of random starts for greedy phase
#' @return Comprehensive Result 4 analysis results
run_result4_for_tree <- function(dist_obj, subset_size, maximize = TRUE, n_greedy_starts = 1) {
  
  tree_name <- ifelse("tree_name" %in% names(dist_obj), 
                     dist_obj$tree_name, 
                     paste0("tree_", length(dist_obj$tip_labels), "tips"))
  
  n_tips <- length(dist_obj$tip_labels)
  total_combinations <- choose(n_tips, subset_size)
  
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat("Result 4 analysis for", tree_name, "\n")
  cat("Tree size:", n_tips, "tips\n")
  cat("Subset size:", subset_size, "\n")
  cat("Total combinations:", total_combinations, "\n")
  cat("Optimization:", ifelse(maximize, "Maximization", "Minimization"), "\n")
  cat(paste0(strrep("=", 70), "\n"))
  
  # Step 1: Run heuristic algorithm
  cat("Step 1: Running heuristic algorithm...\n")
  heuristic_result <- run_complete_algorithm(dist_obj, subset_size, maximize, 
                                            n_greedy_starts = n_greedy_starts)
  
  cat("  Heuristic subset:", heuristic_result$final_subset_names, "\n")
  cat("  Heuristic metrics: MinPD =", heuristic_result$final_metrics$MinPD,
      "MeanPD =", heuristic_result$final_metrics$MeanPD,
      "MeanNND =", heuristic_result$final_metrics$MeanNND, "\n")
  
  # Step 2: Find exact optimum
  cat("\nStep 2: Finding exact optimum...\n")
  exact_result <- find_exact_optimum(dist_obj, subset_size, maximize)
  
  cat("  Exact optimum subset:", exact_result$subset_names, "\n")
  cat("  Exact optimum metrics: MinPD =", exact_result$metrics$MinPD,
      "MeanPD =", exact_result$metrics$MeanPD,
      "MeanNND =", exact_result$metrics$MeanNND, "\n")
  
  # Step 3: Compare heuristic with exact optimum
  cat("\nStep 3: Comparing heuristic with exact optimum...\n")
  comparison <- compare_heuristic_vs_exact(dist_obj, heuristic_result, subset_size, maximize)
  
  # Step 4: Evaluate heuristic performance
  cat("\nStep 4: Evaluating heuristic performance...\n")
  evaluation <- evaluate_heuristic_performance(dist_obj, subset_size, maximize, n_greedy_starts)
  
  # Step 5: Create summary
  summary <- list(
    tree_name = tree_name,
    tree_n_tips = n_tips,
    subset_size = subset_size,
    total_combinations = total_combinations,
    maximize = maximize,
    heuristic_result = heuristic_result,
    exact_result = exact_result,
    comparison = comparison,
    evaluation = evaluation,
    is_exact_match = comparison$is_exact_match
  )
  
  # Step 6: Print results
  cat("\nStep 5: Results summary\n")
  cat(paste0(strrep("-", 60), "\n"))
  
  # Print comparison summary
  cat("Comparison of heuristic vs exact optimum:\n\n")
  print(comparison$comparison_summary)
  
  cat("\nHeuristic found exact optimum?", comparison$is_exact_match, "\n")
  
  if (!comparison$is_exact_match) {
    cat("\nGaps between heuristic and exact optimum:\n")
    cat("  MinPD gap:  ", sprintf("%+.3f", comparison$gaps$MinPD), 
        "(", sprintf("%+.1f%%", comparison$rel_gaps$MinPD), ")\n")
    cat("  MeanPD gap: ", sprintf("%+.3f", comparison$gaps$MeanPD), 
        "(", sprintf("%+.1f%%", comparison$rel_gaps$MeanPD), ")\n")
    cat("  MeanNND gap:", sprintf("%+.3f", comparison$gaps$MeanNND), 
        "(", sprintf("%+.1f%%", comparison$rel_gaps$MeanNND), ")\n")
  }
  
  # Print rank analysis
  if (!is.null(evaluation$rank_analysis)) {
    cat("\nRank analysis:\n")
    cat("  Better subsets found:", evaluation$rank_analysis$better_count, "\n")
    cat("  Equal subsets found: ", evaluation$rank_analysis$equal_count, "\n")
    cat("  Evaluated subsets:   ", evaluation$rank_analysis$evaluated_count, "\n")
    cat("  Rank of heuristic:   ", evaluation$rank_analysis$rank, "\n")
    cat("  Percentile:          ", round(evaluation$rank_analysis$percentile, 2), "%\n")
  }
  
  cat(paste0(strrep("-", 60), "\n"))
  # Step 7: Conclusions
  cat("\nStep 6: Conclusions\n")
  cat(paste0(strrep("-", 60), "\n"))
  
  if (comparison$is_exact_match) {
    cat("The heuristic algorithm found the EXACT OPTIMUM!\n")
    cat("This provides strong validation of the heuristic approach.\n")
  } else {
    cat("The heuristic algorithm did NOT find the exact optimum.\n")
    
    # Check how close it is
    max_rel_gap <- max(abs(unlist(comparison$rel_gaps)), na.rm = TRUE)
    
    if (max_rel_gap < 1) {
      cat("However, the heuristic is VERY CLOSE to the optimum (max gap < 1%).\n")
      cat("This suggests the heuristic is highly effective.\n")
    } else if (max_rel_gap < 5) {
      cat("The heuristic is CLOSE to the optimum (max gap < 5%).\n")
      cat("This suggests the heuristic is effective.\n")
    } else if (max_rel_gap < 10) {
      cat("The heuristic is MODERATELY close to the optimum (max gap < 10%).\n")
      cat("This suggests the heuristic is reasonably effective.\n")
    } else {
      cat("The heuristic is NOT very close to the optimum (max gap >= 10%).\n")
      cat("This suggests the heuristic may need improvement.\n")
    }
    
    # Check rank
    if (!is.null(evaluation$rank_analysis)) {
      percentile <- evaluation$rank_analysis$percentile
      
      if (percentile >= 99) {
        cat("The heuristic is in the TOP 1% of all possible subsets.\n")
        cat("This is excellent performance.\n")
      } else if (percentile >= 95) {
        cat("The heuristic is in the TOP 5% of all possible subsets.\n")
        cat("This is very good performance.\n")
      } else if (percentile >= 90) {
        cat("The heuristic is in the TOP 10% of all possible subsets.\n")
        cat("This is good performance.\n")
      } else if (percentile >= 75) {
        cat("The heuristic is in the TOP 25% of all possible subsets.\n")
        cat("This is acceptable performance.\n")
      } else {
        cat("The heuristic is not in the top quartile of all possible subsets.\n")
        cat("This suggests room for improvement.\n")
      }
    }
  }
  
  return(summary)
}

#' Run Result 4 analysis for small trees (balanced and ladder)
#'
#' @param dist_objs List of distance objects
#' @param cfg Configuration list
#' @return Comprehensive Result 4 results
run_result4_analysis <- function(dist_objs, cfg) {
  cat(paste0("\n", strrep("=", 70), "\n"))
  cat("RESULT 4: Comparing heuristic with exact optimum\n")
  cat(paste0(strrep("=", 70), "\n"))
  
  results <- list()
  summary_dfs <- list()
  
  # Process small balanced tree
  if ("small_balanced" %in% names(dist_objs)) {
    dist_obj <- dist_objs$small_balanced
    
    cat("\n>>> Processing small balanced tree (n =", cfg$small_n, ")\n")
    cat("Subset size:", cfg$subset_small, "\n")
    cat("Total combinations: choose(", cfg$small_n, ",", cfg$subset_small, ") =", 
        choose(cfg$small_n, cfg$subset_small), "\n")
    
    # Check if exhaustive search is feasible
    total_combinations <- choose(cfg$small_n, cfg$subset_small)
    if (total_combinations > 1000000) {
      cat("Warning: Total combinations (", total_combinations, ") is very large.\n")
      cat("Exhaustive search may be computationally expensive.\n")
      cat("Consider using a smaller tree or subset size.\n")
    }
    
    # Run analysis
    result <- run_result4_for_tree(
      dist_obj,
      subset_size = cfg$subset_small,
      maximize = TRUE,
      n_greedy_starts = 1
    )
    
    results$small_balanced <- result
    
    # Create summary
    summary_df <- data.frame(
      Tree = result$tree_name,
      N_Tips = result$tree_n_tips,
      Subset_Size = result$subset_size,
      Total_Combinations = result$total_combinations,
      Exact_Match = result$is_exact_match,
      Heuristic_MinPD = result$heuristic_result$final_metrics$MinPD,
      Exact_MinPD = result$exact_result$metrics$MinPD,
      MinPD_Gap = result$comparison$gaps$MinPD,
      MinPD_Rel_Gap = result$comparison$rel_gaps$MinPD,
      Heuristic_MeanPD = result$heuristic_result$final_metrics$MeanPD,
      Exact_MeanPD = result$exact_result$metrics$MeanPD,
      MeanPD_Gap = result$comparison$gaps$MeanPD,
      MeanPD_Rel_Gap = result$comparison$rel_gaps$MeanPD,
      Heuristic_MeanNND = result$heuristic_result$final_metrics$MeanNND,
      Exact_MeanNND = result$exact_result$metrics$MeanNND,
      MeanNND_Gap = result$comparison$gaps$MeanNND,
      MeanNND_Rel_Gap = result$comparison$rel_gaps$MeanNND,
      stringsAsFactors = FALSE
    )
    
    # Add rank information if available
    if (!is.null(result$evaluation$rank_analysis)) {
      summary_df$Rank <- result$evaluation$rank_analysis$rank
      summary_df$Percentile <- result$evaluation$rank_analysis$percentile
      summary_df$Evaluated_Subsets <- result$evaluation$rank_analysis$evaluated_count
    }
    
    summary_dfs$small_balanced <- summary_df
  }
  
  # Process small ladder tree
  if ("small_ladder" %in% names(dist_objs)) {
    dist_obj <- dist_objs$small_ladder
    
    cat("\n>>> Processing small ladder tree (n =", cfg$small_n, ")\n")
    cat("Subset size:", cfg$subset_small, "\n")
    cat("Total combinations: choose(", cfg$small_n, ",", cfg$subset_small, ") =", 
        choose(cfg$small_n, cfg$subset_small), "\n")
    
    # Check if exhaustive search is feasible
    total_combinations <- choose(cfg$small_n, cfg$subset_small)
    if (total_combinations > 1000000) {
      cat("Warning: Total combinations (", total_combinations, ") is very large.\n")
      cat("Exhaustive search may be computationally expensive.\n")
      cat("Consider using a smaller tree or subset size.\n")
    }
    
    # Run analysis
    result <- run_result4_for_tree(
      dist_obj,
      subset_size = cfg$subset_small,
      maximize = TRUE,
      n_greedy_starts = 1
    )
    
    results$small_ladder <- result
    
    # Create summary
    summary_df <- data.frame(
      Tree = result$tree_name,
      N_Tips = result$tree_n_tips,
      Subset_Size = result$subset_size,
      Total_Combinations = result$total_combinations,
      Exact_Match = result$is_exact_match,
      Heuristic_MinPD = result$heuristic_result$final_metrics$MinPD,
      Exact_MinPD = result$exact_result$metrics$MinPD,
      MinPD_Gap = result$comparison$gaps$MinPD,
      MinPD_Rel_Gap = result$comparison$rel_gaps$MinPD,
      Heuristic_MeanPD = result$heuristic_result$final_metrics$MeanPD,
      Exact_MeanPD = result$exact_result$metrics$MeanPD,
      MeanPD_Gap = result$comparison$gaps$MeanPD,
      MeanPD_Rel_Gap = result$comparison$rel_gaps$MeanPD,
      Heuristic_MeanNND = result$heuristic_result$final_metrics$MeanNND,
      Exact_MeanNND = result$exact_result$metrics$MeanNND,
      MeanNND_Gap = result$comparison$gaps$MeanNND,
      MeanNND_Rel_Gap = result$comparison$rel_gaps$MeanNND,
      stringsAsFactors = FALSE
    )
    
    # Add rank information if available
    if (!is.null(result$evaluation$rank_analysis)) {
      summary_df$Rank <- result$evaluation$rank_analysis$rank
      summary_df$Percentile <- result$evaluation$rank_analysis$percentile
      summary_df$Evaluated_Subsets <- result$evaluation$rank_analysis$evaluated_count
    }
    
    summary_dfs$small_ladder <- summary_df
  }
  
  # Combine all summaries
  if (length(summary_dfs) > 0) {
    combined_summary <- do.call(rbind, summary_dfs)
    rownames(combined_summary) <- NULL
    
    cat(paste0("\n", strrep("=", 70), "\n"))
    cat("OVERALL SUMMARY\n")
    cat(paste0(strrep("=", 70), "\n"))
    print(combined_summary)
    
    return(list(
      results = results,
      summary = combined_summary
    ))
  } else {
    cat("No small trees found for Result 4 analysis.\n")
    cat("Result 4 analysis requires small trees for exhaustive search.\n")
    return(NULL)
  }
}

#' Save Result 4 results to files
#'
#' @param result4_results Results from run_result4_analysis
#' @param cfg Configuration list
save_result4 <- function(result4_results, cfg) {
  cat("\nSaving Result 4 results...\n")
  
  # Create output directory if it doesn't exist
  output_dir <- cfg$tables_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save summary
  if (!is.null(result4_results$summary)) {
    summary_file <- file.path(output_dir, "result4_summary.csv")
    write.csv(result4_results$summary, summary_file, row.names = FALSE)
    cat("  Saved summary to:", summary_file, "\n")
  }
  
  # Save detailed comparison for each tree
  if (!is.null(result4_results$results)) {
    for (tree_name in names(result4_results$results)) {
      result <- result4_results$results[[tree_name]]
      
      # Save heuristic result
      heuristic_file <- file.path(output_dir, paste0("result4_heuristic_", tree_name, ".csv"))
      heuristic_df <- data.frame(
        Tree = result$tree_name,
        Subset_Size = result$subset_size,
        Tip_Index = result$heuristic_result$final_subset,
        Tip_Name = result$heuristic_result$final_subset_names,
        MinPD = result$heuristic_result$final_metrics$MinPD,
        MeanPD = result$heuristic_result$final_metrics$MeanPD,
        MeanNND = result$heuristic_result$final_metrics$MeanNND,
        stringsAsFactors = FALSE
      )
      write.csv(heuristic_df, heuristic_file, row.names = FALSE)
      
      # Save exact result
      exact_file <- file.path(output_dir, paste0("result4_exact_", tree_name, ".csv"))
      exact_df <- data.frame(
        Tree = result$tree_name,
        Subset_Size = result$subset_size,
        Tip_Index = result$exact_result$subset,
        Tip_Name = result$exact_result$subset_names,
        MinPD = result$exact_result$metrics$MinPD,
        MeanPD = result$exact_result$metrics$MeanPD,
        MeanNND = result$exact_result$metrics$MeanNND,
        Combination_Index = result$exact_result$combination_index,
        Total_Combinations = result$exact_result$total_combinations,
        stringsAsFactors = FALSE
      )
      write.csv(exact_df, exact_file, row.names = FALSE)
      
      # Save comparison
      comparison_file <- file.path(output_dir, paste0("result4_comparison_", tree_name, ".csv"))
      write.csv(result$comparison$comparison_summary, comparison_file, row.names = FALSE)
      
      # Save rank analysis if available
      if (!is.null(result$evaluation$rank_analysis)) {
        rank_file <- file.path(output_dir, paste0("result4_rank_analysis_", tree_name, ".csv"))
        rank_df <- data.frame(
          Tree = result$tree_name,
          Better_Subsets = result$evaluation$rank_analysis$better_count,
          Equal_Subsets = result$evaluation$rank_analysis$equal_count,
          Evaluated_Subsets = result$evaluation$rank_analysis$evaluated_count,
          Rank = result$evaluation$rank_analysis$rank,
          Percentile = result$evaluation$rank_analysis$percentile,
          stringsAsFactors = FALSE
        )
        write.csv(rank_df, rank_file, row.names = FALSE)
      }
    }
  }
  
  cat("Result 4 results saved to:", output_dir, "\n")
}

#' Test Result 4 analysis
test_result4 <- function() {
  # Load required libraries and functions
  library(ape)
  source("config/analysis_config.R")
  source("distance_metrics.R")
  source("objective_compare.R")
  source("subset_greedy.R")
  source("subset_exchange.R")
  source("subset_exhaustive.R")
  
  # Create a test configuration
  test_cfg <- list(
    small_n = 15,  # Small for testing
    subset_small = 5
  )
  
  # Create a small test tree
  cat("Creating test tree...\n")
  test_tree <- rtree(test_cfg$small_n)
  test_tree$tip.label <- paste0("sp", 1:test_cfg$small_n)
  dist_obj <- create_distance_object(test_tree)
  dist_obj$tree_name <- "test_small_tree"
  
  # Create distance objects list
  dist_objs <- list(small_balanced = dist_obj)
  
  # Run Result 4 analysis
  cat("\nRunning Result 4 analysis...\n")
  result4 <- run_result4_analysis(dist_objs, test_cfg)
  
  # Save results
  test_cfg$tables_dir <- "test_output"
  if (!dir.exists(test_cfg$tables_dir)) {
    dir.create(test_cfg$tables_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  save_result4(result4, test_cfg)
  
  return(result4)
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_result <- test_result4()
  cat("\nTest completed successfully.\n")
}
