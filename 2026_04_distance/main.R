# Main program for phylogenetic dispersed subset analysis
# This script orchestrates the entire analysis pipeline

# Clear workspace
rm(list = ls())

# Set working directory to project root
setwd("/home/huangr/projects/2026_04_distance")

# Load configuration and utility functions
source("config/analysis_config.R")
source("R/utils_io.R")

#' Main analysis function
#'
#' @param cfg Configuration list (optional)
#' @return Comprehensive results from all analyses
run_analysis <- function(cfg = NULL) {
  cat(paste0("\n", strrep("=", 80), "\n"))
  cat("PHYLOGENETIC DISPERSED SUBSET ANALYSIS\n")
  cat(paste0(strrep("=", 80), "\n"))
  
  # Use default configuration if not provided
  if (is.null(cfg)) {
    # Load configuration from file
    source("config/analysis_config.R")
    # cfg is now available in the global environment after sourcing
    # We need to get it from the parent frame
    cfg <- get("cfg", envir = .GlobalEnv)
  }
  
  # Set random seed for reproducibility
  set.seed(cfg$seed)
  
  # Create output directories
  cat("\n1. Setting up output directories...\n")
  dirs <- create_output_dirs()
  cfg$tables_dir <- dirs$tables
  cfg$figures_dir <- dirs$figures
  
  # Create log file
  log_con <- create_log()
  write_log(log_con, "Analysis started")
  write_log(log_con, paste("Configuration:", paste(names(cfg), collapse = ", ")))
  
  # Load required packages
  write_log(log_con, "Loading required packages...")
  load_required_packages(install_missing = TRUE)
  
  # Load all analysis modules
  write_log(log_con, "Loading analysis modules...")
  source_modules()
  
  # Generate trees
  write_log(log_con, "Generating phylogenetic trees...")
  cat("\n2. Generating phylogenetic trees...\n")
  trees <- generate_all_trees(cfg)
  
  # Create distance objects
  write_log(log_con, "Creating distance objects...")
  cat("\n3. Creating distance objects...\n")
  dist_objs <- create_distance_objects(trees)
  
  # Save distance objects for future use
  write_log(log_con, "Saving distance objects...")
  save_rds(dist_objs, "distance_objects", dirs$rds)
  
  # Initialize results container
  all_results <- list()
  
  # Run Result 1: Comparison with random null distribution
  if (cfg$run_result1) {
    write_log(log_con, "Starting Result 1 analysis...")
    cat(paste0("\n", strrep("-", 70), "\n"))
    cat("RESULT 1: Comparison with Random Null Distribution\n")
    cat(paste0(strrep("-", 70), "\n"))
    
    result1 <- run_result1_analysis(dist_objs, cfg)
    all_results$result1 <- result1
    
    # Save results
    save_result1(result1, cfg)
    save_rds(result1, "result1", dirs$rds)
    
    # Create figures
    write_log(log_con, "Creating Result 1 figures...")
    for (tree_name in names(result1$results)) {
      create_result1_figures(result1, tree_name, dirs$figures)
    }
    
    write_log(log_con, "Result 1 analysis completed")
  }
  
  # Run Result 2: Design component analysis
  if (cfg$run_result2) {
    write_log(log_con, "Starting Result 2 analysis...")
    cat(paste0("\n", strrep("-", 70), "\n"))
    cat("RESULT 2: Design Component Analysis\n")
    cat(paste0(strrep("-", 70), "\n"))
    
    result2 <- run_result2_analysis(dist_objs, cfg)
    all_results$result2 <- result2
    
    # Save results
    save_result2(result2, cfg)
    save_rds(result2, "result2", dirs$rds)
    
    # Create figures
    write_log(log_con, "Creating Result 2 figures...")
    for (tree_name in names(result2$results)) {
      create_result2_figures(result2, tree_name, dirs$figures)
    }
    
    write_log(log_con, "Result 2 analysis completed")
  }
  
  # Run Result 4: Heuristic vs exact optimum (only for small tree)
  if (cfg$run_result4) {
    write_log(log_con, "Starting Result 4 analysis...")
    cat(paste0("\n", strrep("-", 70), "\n"))
    cat("RESULT 4: Heuristic vs Exact Optimum\n")
    cat(paste0(strrep("-", 70), "\n"))
    
    result4 <- run_result4_analysis(dist_objs, cfg)
    all_results$result4 <- result4
    
    # Save results
    save_result4(result4, cfg)
    save_rds(result4, "result4", dirs$rds)
    
    # Create figures
    write_log(log_con, "Creating Result 4 figures...")
    create_result4_figures(result4, dirs$figures)
    
    write_log(log_con, "Result 4 analysis completed")
  }
  
  # Run Result 3: Trait-level consequences (runs last as per design)
  if (cfg$run_result3) {
    write_log(log_con, "Starting Result 3 analysis...")
    cat(paste0("\n", strrep("-", 70), "\n"))
    cat("RESULT 3: Trait-level Consequences\n")
    cat(paste0(strrep("-", 70), "\n"))
    
    result3 <- run_result3_analysis(dist_objs, cfg)
    all_results$result3 <- result3
    
    # Save results
    save_result3(result3, cfg)
    save_rds(result3, "result3", dirs$rds)
    
    # Create figures
    write_log(log_con, "Creating Result 3 figures...")
    for (tree_name in names(result3$results)) {
      create_result3_figures(result3, tree_name, dirs$figures)
    }
    
    write_log(log_con, "Result 3 analysis completed")
  }
  
  # Create summary figures
  write_log(log_con, "Creating summary figures...")
  cat(paste0("\n", strrep("-", 70), "\n"))
  cat("CREATING SUMMARY FIGURES\n")
  cat(paste0(strrep("-", 70), "\n"))
  
  summary_plots <- create_summary_figures(all_results, cfg, dirs$figures)
  save_rds(summary_plots, "summary_plots", dirs$rds)
  
  # Create summary report
  write_log(log_con, "Creating summary report...")
  report_file <- create_summary_report(all_results, cfg)
  
  # Save all results
  write_log(log_con, "Saving all results...")
  save_rds(all_results, "all_results", dirs$rds)
  
  # Log completion
  write_log(log_con, "Analysis completed successfully")
  close_log(log_con)
  
  # Print completion message
  cat(paste0("\n", strrep("=", 80), "\n"))
  cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
  cat(paste0(strrep("=", 80), "\n"))
  
  cat("\nOutput directories:\n")
  cat("  Tables:", dirs$tables, "\n")
  cat("  Figures:", dirs$figures, "\n")
  cat("  RDS objects:", dirs$rds, "\n")
  cat("  Logs:", dirs$logs, "\n")
  
  cat("\nSummary report:", report_file, "\n")
  
  return(all_results)
}

#' Source all analysis modules
#'
#' @return NULL
source_modules <- function() {
  cat("Loading analysis modules...\n")
  
  # Core modules
  source("R/tree_generators.R")
  source("R/distance_metrics.R")
  source("R/objective_compare.R")
  source("R/subset_greedy.R")
  source("R/subset_exchange.R")
  source("R/subset_random.R")
  source("R/subset_exhaustive.R")
  source("R/trait_simulation.R")
  source("R/signal_metrics.R")
  
  # Analysis modules
  source("R/result1_analysis.R")
  source("R/result2_analysis.R")
  source("R/result3_analysis.R")
  source("R/result4_analysis.R")
  
  # Visualization module
  source("R/plotting.R")
  
  cat("All modules loaded successfully.\n")
}

#' Create distance objects from trees
#'
#' @param trees List of phylogenetic trees
#' @return List of distance objects
create_distance_objects <- function(trees) {
  dist_objs <- list()
  
  for (tree_name in names(trees)) {
    tree <- trees[[tree_name]]
    
    dist_obj <- create_distance_object(tree)
    dist_obj$tree_name <- tree_name
    
    dist_objs[[tree_name]] <- dist_obj
  }
  
  return(dist_objs)
}

#' Run a quick test of the analysis pipeline
#'
#' @return Test results
test_analysis <- function() {
  cat(paste0("\n", strrep("=", 80), "\n"))
  cat("TESTING ANALYSIS PIPELINE\n")
  cat(paste0(strrep("=", 80), "\n"))
  
  # Create test configuration
  source("config/analysis_config.R")
  test_cfg <- cfg
  
  # Modify for testing (smaller sizes for speed)
  test_cfg$large_n <- 64
  test_cfg$small_n <- 16
  test_cfg$subset_large <- 10
  test_cfg$subset_small <- 5
  test_cfg$null_reps_large <- 100
  test_cfg$null_reps_small <- 100
  test_cfg$trait_reps <- 10
  test_cfg$random_subset_reps_for_trait <- 20
  test_cfg$exhaustive_small <- FALSE  # Don't run exhaustive search in test
  
  cat("\nTest configuration:\n")
  print(test_cfg)
  
  # Run analysis
  cat("\nRunning test analysis...\n")
  test_results <- run_analysis(test_cfg)
  
  cat(paste0("\n", strrep("=", 80), "\n"))
  cat("TEST COMPLETED SUCCESSFULLY\n")
  cat(paste0(strrep("=", 80), "\n"))
  
  return(test_results)
}

#' Run individual result analysis
#'
#' @param result_number Which result to run (1, 2, 3, or 4)
#' @param cfg Configuration list (optional)
#' @return Results for the specified analysis
run_single_result <- function(result_number, cfg = NULL) {
  if (is.null(cfg)) {
    source("config/analysis_config.R")
    cfg <- get("cfg", envir = .GlobalEnv)
  }
  
  set.seed(cfg$seed)
  
  # Create output directories
  dirs <- create_output_dirs()
  cfg$tables_dir <- dirs$tables
  cfg$figures_dir <- dirs$figures
  
  # Load required packages
  load_required_packages(install_missing = TRUE)
  
  # Load all analysis modules
  source_modules()
  
  # Generate trees and distance objects
  trees <- generate_all_trees(cfg)
  dist_objs <- create_distance_objects(trees)
  
  # Run specified result
  if (result_number == 1) {
    cat("\nRunning Result 1 only...\n")
    result <- run_result1_analysis(dist_objs, cfg)
    save_result1(result, cfg)
    
    # Create figures
    for (tree_name in names(result$results)) {
      create_result1_figures(result, tree_name, dirs$figures)
    }
    
  } else if (result_number == 2) {
    cat("\nRunning Result 2 only...\n")
    result <- run_result2_analysis(dist_objs, cfg)
    save_result2(result, cfg)
    
    # Create figures
    for (tree_name in names(result$results)) {
      create_result2_figures(result, tree_name, dirs$figures)
    }
    
  } else if (result_number == 3) {
    cat("\nRunning Result 3 only...\n")
    result <- run_result3_analysis(dist_objs, cfg)
    save_result3(result, cfg)
    
    # Create figures
    for (tree_name in names(result$results)) {
      create_result3_figures(result, tree_name, dirs$figures)
    }
    
  } else if (result_number == 4) {
    cat("\nRunning Result 4 only...\n")
    result <- run_result4_analysis(dist_objs, cfg)
    save_result4(result, cfg)
    
    # Create figures
    create_result4_figures(result, dirs$figures)
    
  } else {
    stop("Invalid result number. Must be 1, 2, 3, or 4.")
  }
  
  return(result)
}

# Main execution
if (sys.nframe() == 0) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    # Run full analysis by default
    cat("No arguments provided. Running full analysis...\n")
    results <- run_analysis()
    
  } else if (args[1] == "test") {
    # Run test
    test_results <- test_analysis()
    
  } else if (args[1] == "result1") {
    # Run Result 1 only
    result1 <- run_single_result(1)
    
  } else if (args[1] == "result2") {
    # Run Result 2 only
    result2 <- run_single_result(2)
    
  } else if (args[1] == "result3") {
    # Run Result 3 only
    result3 <- run_single_result(3)
    
  } else if (args[1] == "result4") {
    # Run Result 4 only
    result4 <- run_single_result(4)
    
  } else if (args[1] == "help") {
    # Show help
    cat("\nUsage: Rscript main.R [option]\n")
    cat("\nOptions:\n")
    cat("  (no argument)    Run full analysis\n")
    cat("  test             Run test analysis (smaller sizes)\n")
    cat("  result1          Run Result 1 only\n")
    cat("  result2          Run Result 2 only\n")
    cat("  result3          Run Result 3 only\n")
    cat("  result4          Run Result 4 only\n")
    cat("  help             Show this help message\n")
    
  } else {
    cat("Unknown option:", args[1], "\n")
    cat("Use 'help' to see available options.\n")
  }
}
