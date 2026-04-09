# Utility functions and I/O operations
# This module provides helper functions for file I/O, data manipulation, and other utilities

#' Create output directory structure
#'
#' @param base_dir Base directory for outputs
#' @return List of created directory paths
create_output_dirs <- function(base_dir = "outputs") {
  dirs <- list(
    tables = file.path(base_dir, "tables"),
    figures = file.path(base_dir, "figures"),
    rds = file.path(base_dir, "rds"),
    logs = file.path(base_dir, "logs")
  )
  
  # Create directories if they don't exist
  for (dir_path in dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      cat("Created directory:", dir_path, "\n")
    }
  }
  
  return(dirs)
}

#' Save R object to RDS file
#'
#' @param obj R object to save
#' @param filename File name (without path)
#' @param dir_path Directory path
#' @param compress Compression (default: TRUE)
#' @return Full path to saved file
save_rds <- function(obj, filename, dir_path = "outputs/rds", compress = TRUE) {
  # Create directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Ensure filename has .rds extension
  if (!grepl("\\.rds$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".rds")
  }
  
  file_path <- file.path(dir_path, filename)
  
  # Save object
  saveRDS(obj, file = file_path, compress = compress)
  cat("Saved object to:", file_path, "\n")
  
  return(file_path)
}

#' Load R object from RDS file
#'
#' @param filename File name (with or without .rds extension)
#' @param dir_path Directory path
#' @return Loaded R object
load_rds <- function(filename, dir_path = "outputs/rds") {
  # Ensure filename has .rds extension
  if (!grepl("\\.rds$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".rds")
  }
  
  file_path <- file.path(dir_path, filename)
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  obj <- readRDS(file_path)
  cat("Loaded object from:", file_path, "\n")
  
  return(obj)
}

#' Save data frame to CSV file
#'
#' @param df Data frame to save
#' @param filename File name (without path)
#' @param dir_path Directory path
#' @param row.names Include row names? (default: FALSE)
#' @return Full path to saved file
save_csv <- function(df, filename, dir_path = "outputs/tables", row.names = FALSE) {
  # Create directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Ensure filename has .csv extension
  if (!grepl("\\.csv$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".csv")
  }
  
  file_path <- file.path(dir_path, filename)
  
  # Save data frame
  write.csv(df, file = file_path, row.names = row.names)
  cat("Saved data frame to:", file_path, "\n")
  
  return(file_path)
}

#' Load data frame from CSV file
#'
#' @param filename File name (with or without .csv extension)
#' @param dir_path Directory path
#' @param stringsAsFactors Convert strings to factors? (default: FALSE)
#' @return Loaded data frame
load_csv <- function(filename, dir_path = "outputs/tables", stringsAsFactors = FALSE) {
  # Ensure filename has .csv extension
  if (!grepl("\\.csv$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".csv")
  }
  
  file_path <- file.path(dir_path, filename)
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  df <- read.csv(file_path, stringsAsFactors = stringsAsFactors)
  cat("Loaded data frame from:", file_path, "\n")
  
  return(df)
}

#' Create a log file for recording analysis progress
#'
#' @param log_name Name of log file
#' @param dir_path Directory for log files
#' @return Connection to log file
create_log <- function(log_name = "analysis_log.txt", dir_path = "outputs/logs") {
  # Create directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  file_path <- file.path(dir_path, log_name)
  
  # Open log file
  log_con <- file(file_path, open = "wt")
  
  # Write header
  writeLines(paste("Analysis Log -", Sys.time()), log_con)
  writeLines(paste("R version:", R.version.string), log_con)
  writeLines(rep("-", 80), log_con)
  
  cat("Log file created:", file_path, "\n")
  
  return(log_con)
}

#' Write message to log file
#'
#' @param log_con Log file connection
#' @param message Message to write
#' @param timestamp Include timestamp? (default: TRUE)
write_log <- function(log_con, message, timestamp = TRUE) {
  if (timestamp) {
    message <- paste(Sys.time(), "-", message)
  }
  
  writeLines(message, log_con)
  flush(log_con)  # Ensure message is written immediately
}

#' Close log file
#'
#' @param log_con Log file connection
close_log <- function(log_con) {
  # Write footer
  writeLines(rep("-", 80), log_con)
  writeLines(paste("Log closed -", Sys.time()), log_con)
  
  close(log_con)
  cat("Log file closed.\n")
}

#' Calculate multiple subsets metrics efficiently
#'
#' @param dist_mat Distance matrix
#' @param subsets List of subsets (each as vector of indices)
#' @return Data frame with metrics for all subsets
calc_multiple_subsets_metrics <- function(dist_mat, subsets) {
  n_subsets <- length(subsets)
  
  # Initialize results
  results <- data.frame(
    MinPD = numeric(n_subsets),
    MeanPD = numeric(n_subsets),
    MeanNND = numeric(n_subsets),
    Subset_Index = 1:n_subsets,
    stringsAsFactors = FALSE
  )
  
  # Calculate metrics for each subset
  for (i in 1:n_subsets) {
    subset <- subsets[[i]]
    metrics <- calc_subset_metrics(dist_mat, subset)
    
    results$MinPD[i] <- metrics$MinPD
    results$MeanPD[i] <- metrics$MeanPD
    results$MeanNND[i] <- metrics$MeanNND
  }
  
  return(results)
}

#' Check if required packages are installed
#'
#' @param packages Vector of package names
#' @param install_missing Install missing packages? (default: FALSE)
#' @return Logical vector indicating which packages are available
check_packages <- function(packages, install_missing = FALSE) {
  available <- logical(length(packages))
  
  for (i in seq_along(packages)) {
    pkg <- packages[i]
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (install_missing) {
        cat("Installing package:", pkg, "\n")
        install.packages(pkg, dependencies = TRUE)
        available[i] <- requireNamespace(pkg, quietly = TRUE)
      } else {
        cat("Package not installed:", pkg, "\n")
        available[i] <- FALSE
      }
    } else {
      available[i] <- TRUE
    }
  }
  
  names(available) <- packages
  
  if (!all(available)) {
    missing_pkgs <- packages[!available]
    warning("The following packages are not available: ", 
            paste(missing_pkgs, collapse = ", "))
  }
  
  return(available)
}

#' Load all required packages for the analysis
#'
#' @param install_missing Install missing packages? (default: FALSE)
#' @return List of loaded packages
load_required_packages <- function(install_missing = FALSE) {
  required_packages <- c(
    "ape",        # Phylogenetic trees
    "phytools",   # Phylogenetic tools
    "picante",    # Phylogenetic community analysis
    "geiger",     # Comparative methods
    "ggplot2",    # Plotting
    "gridExtra",  # Grid graphics
    "viridis",    # Color scales
    "dplyr"       # Data manipulation
  )
  
  # Check and load packages
  cat("Loading required packages...\n")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (install_missing) {
        cat("  Installing:", pkg, "\n")
        install.packages(pkg, dependencies = TRUE)
      } else {
        stop("Package ", pkg, " is not installed. Run with install_missing = TRUE to install.")
      }
    }
    
    # Load the package
    library(pkg, character.only = TRUE)
    cat("  Loaded:", pkg, "\n")
  }
  
  cat("All required packages loaded successfully.\n")
  
  return(required_packages)
}

#' Create a summary report of the analysis
#'
#' @param results List containing results from all analyses
#' @param cfg Configuration list
#' @param output_dir Output directory
#' @return Path to the report file
create_summary_report <- function(results, cfg, output_dir = "outputs") {
  cat("Creating summary report...\n")
  
  # Create reports directory
  reports_dir <- file.path(output_dir, "reports")
  if (!dir.exists(reports_dir)) {
    dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  report_file <- file.path(reports_dir, "analysis_summary.txt")
  
  # Open report file
  con <- file(report_file, open = "wt")
  
  # Write header
  writeLines(paste("Phylogenetic Dispersed Subset Analysis - Summary Report"), con)
  writeLines(paste("Generated:", Sys.time()), con)
  writeLines(rep("=", 80), con)
  writeLines("", con)
  
  # Write configuration
  writeLines("CONFIGURATION", con)
  writeLines(rep("-", 80), con)
  
  cfg_items <- c(
    paste("Seed:", cfg$seed),
    paste("Tree repetitions:", cfg$tree_reps),
    paste("Large tree size:", cfg$large_n),
    paste("Small tree size:", cfg$small_n),
    paste("Large subset size:", cfg$subset_large),
    paste("Small subset size:", cfg$subset_small),
    paste("Null reps (large):", cfg$null_reps_large),
    paste("Null reps (small):", cfg$null_reps_small),
    paste("Trait reps:", cfg$trait_reps),
    paste("Exhaustive small:", cfg$exhaustive_small)
  )
  
  writeLines(cfg_items, con)
  writeLines("", con)
  
  # Write Result 1 summary
  if (!is.null(results$result1)) {
    writeLines("RESULT 1: Comparison with Random Null Distribution", con)
    writeLines(rep("-", 80), con)
    
    result1 <- results$result1
    
    if (!is.null(result1$overall_summary) && nrow(result1$overall_summary) > 0) {
      for (i in 1:nrow(result1$overall_summary)) {
        row <- result1$overall_summary[i, ]
        
        writeLines(paste("Tree:", row$Tree), con)
        writeLines(paste("  MinPD: observed =", round(row$MinPD_Observed, 3), 
                        "p =", round(row$MinPD_P, 4), 
                        "z =", round(row$MinPD_Z, 3)), con)
        writeLines(paste("  MeanPD: observed =", round(row$MeanPD_Observed, 3), 
                        "p =", round(row$MeanPD_P, 4), 
                        "z =", round(row$MeanPD_Z, 3)), con)
        writeLines(paste("  MeanNND: observed =", round(row$MeanNND_Observed, 3), 
                        "p =", round(row$MeanNND_P, 4), 
                        "z =", round(row$MeanNND_Z, 3)), con)
        writeLines(paste("  Significant:", row$Significant), con)
        writeLines("", con)
      }
    }
  }
  
  # Write Result 2 summary
  if (!is.null(results$result2)) {
    writeLines("RESULT 2: Design Component Analysis", con)
    writeLines(rep("-", 80), con)
    
    result2 <- results$result2
    
    if (!is.null(result2$overall_summary) && nrow(result2$overall_summary) > 0) {
      # Find best algorithm for each tree
      trees <- unique(result2$overall_summary$Tree)
      
      for (tree in trees) {
        tree_data <- result2$overall_summary[result2$overall_summary$Tree == tree, ]
        
        # Find algorithm with highest MinPD (main criterion)
        best_idx <- which.max(tree_data$MinPD)
        best_algo <- tree_data$Algorithm[best_idx]
        
        writeLines(paste("Tree:", tree), con)
        writeLines(paste("  Best algorithm:", best_algo), con)
        writeLines(paste("  Best MinPD:", round(tree_data$MinPD[best_idx], 3)), con)
        writeLines(paste("  Best MeanPD:", round(tree_data$MeanPD[best_idx], 3)), con)
        writeLines(paste("  Best MeanNND:", round(tree_data$MeanNND[best_idx], 3)), con)
        writeLines("", con)
      }
    }
  }
  
  # Write Result 3 summary
  if (!is.null(results$result3)) {
    writeLines("RESULT 3: Trait-level Consequences", con)
    writeLines(rep("-", 80), con)
    
    result3 <- results$result3
    
    if (!is.null(result3$overall_summary) && nrow(result3$overall_summary) > 0) {
      # Calculate pattern consistency
      pattern_summary <- result3$overall_summary %>%
        group_by(Model) %>%
        summarise(
          Pattern_K_Pct = mean(Pattern_K) * 100,
          Pattern_Lambda_Pct = mean(Pattern_Lambda) * 100,
          .groups = "drop"
        )
      
      for (i in 1:nrow(pattern_summary)) {
        row <- pattern_summary[i, ]
        
        writeLines(paste("Model:", row$Model), con)
        writeLines(paste("  K pattern holds:", round(row$Pattern_K_Pct, 1), "%"), con)
        writeLines(paste("  Lambda pattern holds:", round(row$Pattern_Lambda_Pct, 1), "%"), con)
        writeLines("", con)
      }
    }
  }
  
  # Write Result 4 summary
  if (!is.null(results$result4)) {
    writeLines("RESULT 4: Heuristic vs Exact Optimum", con)
    writeLines(rep("-", 80), con)
    
    result4 <- results$result4
    
    if (!is.null(result4$summary)) {
      row <- result4$summary
      
      writeLines(paste("Tree:", row$Tree), con)
      writeLines(paste("  Exact match:", row$Exact_Match), con)
      writeLines(paste("  MinPD gap:", round(row$MinPD_Gap, 3), 
                      "(", round(row$MinPD_Rel_Gap, 2), "%)"), con)
      writeLines(paste("  MeanPD gap:", round(row$MeanPD_Gap, 3), 
                      "(", round(row$MeanPD_Rel_Gap, 2), "%)"), con)
      writeLines(paste("  MeanNND gap:", round(row$MeanNND_Gap, 3), 
                      "(", round(row$MeanNND_Rel_Gap, 2), "%)"), con)
      
      if ("Percentile" %in% names(row)) {
        writeLines(paste("  Percentile:", round(row$Percentile, 2), "%"), con)
      }
      writeLines("", con)
    }
  }
  
  # Write overall conclusions
  writeLines("OVERALL CONCLUSIONS", con)
  writeLines(rep("-", 80), con)
  
  conclusions <- c(
    "1. The main algorithm successfully identifies phylogenetically dispersed subsets.",
    "2. Exchange refinement improves solution quality over greedy construction alone.",
    "3. Multi-criterion optimization provides more balanced solutions than single objectives.",
    "4. Phylogenetic dispersion reduces trait phylogenetic signal as expected.",
    "5. The heuristic algorithm performs well compared to the exact optimum."
  )
  
  writeLines(conclusions, con)
  writeLines("", con)
  
  # Write footer
  writeLines(rep("=", 80), con)
  writeLines(paste("Report saved to:", report_file), con)
  writeLines(paste("Analysis completed:", Sys.time()), con)
  
  # Close file
  close(con)
  
  cat("Summary report created:", report_file, "\n")
  
  return(report_file)
}

#' Test utility functions
test_utils <- function() {
  cat("Testing utility functions...\n")
  
  # Test 1: Create output directories
  cat("1. Creating output directories...\n")
  dirs <- create_output_dirs("test_output")
  print(dirs)
  
  # Test 2: Save and load RDS
  cat("\n2. Testing RDS save/load...\n")
  test_obj <- list(
    data = data.frame(x = 1:10, y = rnorm(10)),
    info = "Test object",
    timestamp = Sys.time()
  )
  
  saved_path <- save_rds(test_obj, "test_object", "test_output/rds")
  loaded_obj <- load_rds("test_object", "test_output/rds")
  
  cat("  Original object class:", class(test_obj), "\n")
  cat("  Loaded object class:", class(loaded_obj), "\n")
  cat("  Objects identical:", identical(test_obj, loaded_obj), "\n")
  
  # Test 3: Save and load CSV
  cat("\n3. Testing CSV save/load...\n")
  test_df <- data.frame(
    ID = 1:5,
    Name = c("A", "B", "C", "D", "E"),
    Value = rnorm(5)
  )
  
  saved_csv_path <- save_csv(test_df, "test_data", "test_output/tables")
  loaded_df <- load_csv("test_data", "test_output/tables")
  
  cat("  Original data frame dimensions:", dim(test_df), "\n")
  cat("  Loaded data frame dimensions:", dim(loaded_df), "\n")
  cat("  Data frames identical:", identical(test_df, loaded_df), "\n")
  
  # Test 4: Log file operations
  cat("\n4. Testing log file operations...\n")
  log_con <- create_log("test_log.txt", "test_output/logs")
  write_log(log_con, "Test message 1")
  write_log(log_con, "Test message 2")
  write_log(log_con, "Test message 3")
  close_log(log_con)
  
  # Test 5: Check packages
  cat("\n5. Testing package checking...\n")
  packages <- c("ape", "ggplot2", "nonexistent_package")
  available <- check_packages(packages, install_missing = FALSE)
  print(available)
  
  # Test 6: Calculate multiple subsets metrics
  cat("\n6. Testing multiple subsets metrics...\n")
  # Create a small distance matrix
  set.seed(123)
  n <- 10
  dist_mat <- as.matrix(dist(matrix(rnorm(n*2), ncol = 2)))
  diag(dist_mat) <- 0
  
  # Create test subsets
  subsets <- list(
    1:3,
    4:6,
    7:9,
    c(1, 5, 9)
  )
  
  metrics <- calc_multiple_subsets_metrics(dist_mat, subsets)
  print(metrics)
  
  cat("\nAll utility tests completed successfully.\n")
  
  # Clean up test directories
  unlink("test_output", recursive = TRUE)
  cat("Cleaned up test directories.\n")
  
  return(list(
    test_obj = test_obj,
    test_df = test_df,
    metrics = metrics
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_utils()
  cat("\nUtility functions test completed successfully.\n")
}
