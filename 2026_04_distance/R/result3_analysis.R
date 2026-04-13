# Result 3 analysis: Trait-level consequences of phylogenetic dispersion
# This module implements the analysis for Result 3

#' Prepare subsets for trait analysis
#'
#' @param dist_obj Distance object
#' @param subset_size Subset size
#' @param n_random_subsets Number of random subsets to generate
#' @return A list containing three types of subsets
prepare_subsets_for_trait_analysis <- function(dist_obj, subset_size, n_random_subsets = 100) {
  cat("Preparing subsets for trait analysis...\n")
  
  subsets <- list()
  subset_names <- c()
  
  # 1. Dispersed subset (main method, maximization)
  cat("  1. Generating dispersed subset...\n")
  dispersed_result <- run_complete_algorithm(dist_obj, subset_size, maximize = TRUE)
  subsets$dispersed <- dispersed_result$final_subset_names
  subset_names <- c(subset_names, "dispersed")
  
  # 2. Clustered subset (reverse optimization, minimization)
  cat("  2. Generating clustered subset...\n")
  clustered_result <- run_complete_algorithm(dist_obj, subset_size, maximize = FALSE)
  subsets$clustered <- clustered_result$final_subset_names
  subset_names <- c(subset_names, "clustered")
  
  # 3. Random subsets
  cat("  3. Generating", n_random_subsets, "random subsets...\n")
  random_subsets <- sample_random_subsets(dist_obj, subset_size, n_random_subsets)
  
  # Convert indices to names
  subsets$random <- lapply(random_subsets, function(indices) {
    dist_obj$tip_labels[indices]
  })
  subset_names <- c(subset_names, rep("random", n_random_subsets))
  
  # Flatten random subsets for easier processing
  all_subsets <- list(
    dispersed = subsets$dispersed,
    clustered = subsets$clustered
  )
  
  for (i in 1:n_random_subsets) {
    all_subsets[[paste0("random_", i)]] <- subsets$random[[i]]
  }
  
  return(list(
    subsets = all_subsets,
    subset_names = subset_names,
    subset_types = c("dispersed", "clustered", rep("random", n_random_subsets)),
    n_random_subsets = n_random_subsets
  ))
}

#' Run trait simulation and signal analysis
#'
#' @param tree Full tree
#' @param subset_data Result from prepare_subsets_for_trait_analysis
#' @param n_reps Number of simulation replicates
#' @param ou_alpha_values Vector of alpha values for OU process
#' @return Comprehensive trait analysis results
run_trait_analysis <- function(tree, subset_data, n_reps = 100, 
                               ou_alpha_values = c(0.2, 1, 5)) {
  cat("Running trait simulation and analysis...\n")
  
  # Prepare subsets in the format needed for simulation
  subset_list <- list()
  subset_name_list <- c()
  
  # Add dispersed and clustered
  subset_list$dispersed <- subset_data$subsets$dispersed
  subset_name_list <- c(subset_name_list, "dispersed")
  
  subset_list$clustered <- subset_data$subsets$clustered
  subset_name_list <- c(subset_name_list, "clustered")
  
  # Add ALL random subsets (not just a sample)
  # This is important for proper empirical p-value calculation
  for (i in 1:subset_data$n_random_subsets) {
    subset_list[[paste0("random_", i)]] <- subset_data$subsets[[paste0("random_", i)]]
    subset_name_list <- c(subset_name_list, paste0("random_", i))
  }
  
  # Simulate traits
  cat("  Simulating traits...\n")
  sim_results <- simulate_traits_for_subsets(
    tree, 
    subset_list, 
    subset_name_list, 
    n_reps = n_reps,
    ou_alpha_values = ou_alpha_values
  )
  
  # Calculate phylogenetic signal metrics
  cat("  Calculating phylogenetic signal metrics...\n")
  signal_results <- analyze_simulated_traits_signal(sim_results, tree, n_reps_to_analyze = n_reps)
  
  # Add subset type information
  signal_results$Subset_Type <- NA
  for (i in 1:nrow(signal_results)) {
    subset_name <- signal_results$Subset[i]
    if (subset_name == "dispersed") {
      signal_results$Subset_Type[i] <- "dispersed"
    } else if (subset_name == "clustered") {
      signal_results$Subset_Type[i] <- "clustered"
    } else {
      signal_results$Subset_Type[i] <- "random"
    }
  }
  
  return(list(
    subset_data = subset_data,
    sim_results = sim_results,
    signal_results = signal_results,
    n_reps = n_reps,
    ou_alpha_values = ou_alpha_values
  ))
}

#' Calculate empirical p-values for trait signal comparisons
#'
#' @param signal_df Data frame with signal metrics
#' @param model_name Model name to analyze
#' @param metric_name Metric name ("K", "Lambda", or "MoransI")
#' @return Data frame with empirical p-values for comparisons
calculate_empirical_pvalues <- function(signal_df, model_name, metric_name) {
  # Filter data for this model and metric
  model_data <- signal_df[signal_df$Model == model_name, ]
  
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  
  # Extract values for each subset type
  dispersed_vals <- model_data[[metric_name]][model_data$Subset_Type == "dispersed"]
  clustered_vals <- model_data[[metric_name]][model_data$Subset_Type == "clustered"]
  random_vals <- model_data[[metric_name]][model_data$Subset_Type == "random"]
  
  # Special handling for Moran's I
  # For Moran's I, larger absolute values (both positive and negative) indicate stronger autocorrelation
  # So we need to compare absolute values for Moran's I
  if (metric_name == "MoransI") {
    # For Moran's I, we want to test if dispersed has weaker autocorrelation than clustered/random
    # This means we compare absolute values: |dispersed| < |clustered| and |dispersed| < |random|
    
    # Convert to absolute values for comparison
    dispersed_vals_abs <- abs(dispersed_vals)
    clustered_vals_abs <- abs(clustered_vals)
    random_vals_abs <- abs(random_vals)
    
    # For dispersed vs clustered: test if |dispersed| < |clustered|
    if (length(dispersed_vals) > 0 && length(clustered_vals) > 0) {
      # Count how many clustered have stronger autocorrelation (larger absolute value)
      count_stronger_clustered <- sum(clustered_vals_abs >= dispersed_vals_abs)
      p_dispersed_vs_clustered <- (count_stronger_clustered + 1) / (length(clustered_vals) + 1)
    } else {
      p_dispersed_vs_clustered <- NA
    }
    
    # For dispersed vs random: test if |dispersed| < |random|
    if (length(dispersed_vals) > 0 && length(random_vals) > 0) {
      # Count how many random have stronger autocorrelation (larger absolute value)
      count_stronger_random <- sum(random_vals_abs >= dispersed_vals_abs)
      p_dispersed_vs_random <- (count_stronger_random + 1) / (length(random_vals) + 1)
    } else {
      p_dispersed_vs_random <- NA
    }
    
    # For random vs clustered: test if |random| < |clustered|
    if (length(random_vals) > 0 && length(clustered_vals) > 0) {
      # Count how many clustered have stronger autocorrelation than random
      total_comparisons <- length(random_vals) * length(clustered_vals)
      count_clustered_stronger <- 0
      for (i in 1:length(random_vals)) {
        count_clustered_stronger <- count_clustered_stronger + sum(clustered_vals_abs >= random_vals_abs[i])
      }
      p_random_vs_clustered <- (count_clustered_stronger + 1) / (total_comparisons + 1)
    } else {
      p_random_vs_clustered <- NA
    }
    
  } else {
    # For K and Lambda: regular comparison (higher values = stronger signal)
    
    # For dispersed vs clustered: test if dispersed < clustered
    if (length(dispersed_vals) > 0 && length(clustered_vals) > 0) {
      # One-sided test: dispersed < clustered
      p_dispersed_vs_clustered <- sum(clustered_vals >= dispersed_vals) / length(clustered_vals)
      # Adjust for small sample size
      p_dispersed_vs_clustered <- (sum(clustered_vals >= dispersed_vals) + 1) / (length(clustered_vals) + 1)
    } else {
      p_dispersed_vs_clustered <- NA
    }
    
    # For dispersed vs random: test if dispersed < random
    if (length(dispersed_vals) > 0 && length(random_vals) > 0) {
      # One-sided test: dispersed < random
      p_dispersed_vs_random <- sum(random_vals >= dispersed_vals) / length(random_vals)
      # Adjust for small sample size
      p_dispersed_vs_random <- (sum(random_vals >= dispersed_vals) + 1) / (length(random_vals) + 1)
    } else {
      p_dispersed_vs_random <- NA
    }
    
    # For random vs clustered: test if random < clustered (less important)
    if (length(random_vals) > 0 && length(clustered_vals) > 0) {
      # Count how many clustered values are greater than or equal to random values
      total_comparisons <- length(random_vals) * length(clustered_vals)
      count_clustered_greater <- 0
      for (r_val in random_vals) {
        count_clustered_greater <- count_clustered_greater + sum(clustered_vals >= r_val)
      }
      p_random_vs_clustered <- (count_clustered_greater + 1) / (total_comparisons + 1)
    } else {
      p_random_vs_clustered <- NA
    }
  }
  
  # Calculate means (original values, not absolute)
  mean_dispersed <- mean(dispersed_vals, na.rm = TRUE)
  mean_clustered <- mean(clustered_vals, na.rm = TRUE)
  mean_random <- mean(random_vals, na.rm = TRUE)
  
  return(data.frame(
    Model = model_name,
    Metric = metric_name,
    Mean_Dispersed = mean_dispersed,
    Mean_Clustered = mean_clustered,
    Mean_Random = mean_random,
    P_Dispersed_vs_Clustered = p_dispersed_vs_clustered,
    P_Dispersed_vs_Random = p_dispersed_vs_random,
    P_Random_vs_Clustered = p_random_vs_clustered,
    N_Dispersed = length(dispersed_vals),
    N_Clustered = length(clustered_vals),
    N_Random = length(random_vals),
    stringsAsFactors = FALSE
  ))
}

#' Analyze trait signal results
#'
#' @param trait_analysis_results Results from run_trait_analysis
#' @return Statistical analysis of trait signals
analyze_trait_signals <- function(trait_analysis_results) {
  cat("Analyzing trait signal results...\n")
  
  signal_df <- trait_analysis_results$signal_results
  
  # Remove rows with NA values for key metrics
  signal_df <- signal_df[complete.cases(signal_df$K, signal_df$Lambda, signal_df$MoransI), ]
  
  if (nrow(signal_df) == 0) {
    warning("No complete cases found in signal results")
    return(NULL)
  }
  
  # Summary statistics by model and subset type
  summary_stats <- data.frame()
  
  models <- unique(signal_df$Model)
  subset_types <- unique(signal_df$Subset_Type)
  
  for (model in models) {
    for (subset_type in subset_types) {
      subset_data <- signal_df[signal_df$Model == model & signal_df$Subset_Type == subset_type, ]
      
      if (nrow(subset_data) > 0) {
        summary_stats <- rbind(summary_stats, data.frame(
          Model = model,
          Subset_Type = subset_type,
          N = nrow(subset_data),
          K_mean = mean(subset_data$K, na.rm = TRUE),
          K_sd = sd(subset_data$K, na.rm = TRUE),
          K_min = min(subset_data$K, na.rm = TRUE),
          K_max = max(subset_data$K, na.rm = TRUE),
          Lambda_mean = mean(subset_data$Lambda, na.rm = TRUE),
          Lambda_sd = sd(subset_data$Lambda, na.rm = TRUE),
          Lambda_min = min(subset_data$Lambda, na.rm = TRUE),
          Lambda_max = max(subset_data$Lambda, na.rm = TRUE),
          MoransI_mean = mean(subset_data$MoransI, na.rm = TRUE),
          MoransI_sd = sd(subset_data$MoransI, na.rm = TRUE),
          MoransI_min = min(subset_data$MoransI, na.rm = TRUE),
          MoransI_max = max(subset_data$MoransI, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Calculate empirical p-values for each model and metric
  empirical_pvalues <- data.frame()
  
  for (model in models) {
    for (metric in c("K", "Lambda", "MoransI")) {
      pval_df <- calculate_empirical_pvalues(signal_df, model, metric)
      if (!is.null(pval_df)) {
        empirical_pvalues <- rbind(empirical_pvalues, pval_df)
      }
    }
  }
  
  # Statistical tests
  statistical_tests <- list()
  
  # For each model, test if clustered > random > dispersed
  for (model in models) {
    model_data <- signal_df[signal_df$Model == model, ]
    
    if (nrow(model_data) > 0) {
      # Kruskal-Wallis test for differences among subset types
      if (length(unique(model_data$Subset_Type)) >= 2) {
        kw_test_K <- kruskal.test(K ~ Subset_Type, data = model_data)
        kw_test_Lambda <- kruskal.test(Lambda ~ Subset_Type, data = model_data)
        kw_test_MoransI <- kruskal.test(MoransI ~ Subset_Type, data = model_data)
        
        statistical_tests[[paste0(model, "_Kruskal_Wallis")]] <- list(
          K = kw_test_K,
          Lambda = kw_test_Lambda,
          MoransI = kw_test_MoransI
        )
      }
      
      # Pairwise comparisons
      pairwise_results <- list()
      subset_pairs <- list(
        c("clustered", "random"),
        c("random", "dispersed"),
        c("clustered", "dispersed")
      )
      
      for (pair in subset_pairs) {
        type1 <- pair[1]
        type2 <- pair[2]
        
        data1_K <- model_data$K[model_data$Subset_Type == type1]
        data2_K <- model_data$K[model_data$Subset_Type == type2]
        
        data1_Lambda <- model_data$Lambda[model_data$Subset_Type == type1]
        data2_Lambda <- model_data$Lambda[model_data$Subset_Type == type2]
        
        data1_MoransI <- model_data$MoransI[model_data$Subset_Type == type1]
        data2_MoransI <- model_data$MoransI[model_data$Subset_Type == type2]
        
        if (length(data1_K) > 1 && length(data2_K) > 1) {
          test_K <- wilcox.test(data1_K, data2_K, alternative = "greater")
          test_Lambda <- wilcox.test(data1_Lambda, data2_Lambda, alternative = "greater")
          test_MoransI <- wilcox.test(data1_MoransI, data2_MoransI, alternative = "greater")
          
          pairwise_results[[paste0(type1, "_vs_", type2)]] <- list(
            K = test_K,
            Lambda = test_Lambda,
            MoransI = test_MoransI
          )
        }
      }
      
      statistical_tests[[paste0(model, "_Pairwise")]] <- pairwise_results
    }
  }
  
  # Check if the expected pattern holds: clustered > random > dispersed
  pattern_analysis <- data.frame()
  
  for (model in models) {
    model_summary <- summary_stats[summary_stats$Model == model, ]
    
    if (nrow(model_summary) >= 3) {
      # Get means for each subset type
      clustered_mean_K <- model_summary$K_mean[model_summary$Subset_Type == "clustered"]
      random_mean_K <- model_summary$K_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_K <- model_summary$K_mean[model_summary$Subset_Type == "dispersed"]
      
      clustered_mean_Lambda <- model_summary$Lambda_mean[model_summary$Subset_Type == "clustered"]
      random_mean_Lambda <- model_summary$Lambda_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_Lambda <- model_summary$Lambda_mean[model_summary$Subset_Type == "dispersed"]
      
      clustered_mean_MoransI <- model_summary$MoransI_mean[model_summary$Subset_Type == "clustered"]
      random_mean_MoransI <- model_summary$MoransI_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_MoransI <- model_summary$MoransI_mean[model_summary$Subset_Type == "dispersed"]
      
      # Check pattern
      pattern_K <- (clustered_mean_K > random_mean_K) && (random_mean_K > dispersed_mean_K)
      pattern_Lambda <- (clustered_mean_Lambda > random_mean_Lambda) && (random_mean_Lambda > dispersed_mean_Lambda)
      pattern_MoransI <- (clustered_mean_MoransI > random_mean_MoransI) && (random_mean_MoransI > dispersed_mean_MoransI)
      
      pattern_analysis <- rbind(pattern_analysis, data.frame(
        Model = model,
        Pattern_K = pattern_K,
        Pattern_Lambda = pattern_Lambda,
        Pattern_MoransI = pattern_MoransI,
        Clustered_Mean_K = clustered_mean_K,
        Random_Mean_K = random_mean_K,
        Dispersed_Mean_K = dispersed_mean_K,
        Clustered_Mean_Lambda = clustered_mean_Lambda,
        Random_Mean_Lambda = random_mean_Lambda,
        Dispersed_Mean_Lambda = dispersed_mean_Lambda,
        Clustered_Mean_MoransI = clustered_mean_MoransI,
        Random_Mean_MoransI = random_mean_MoransI,
        Dispersed_Mean_MoransI = dispersed_mean_MoransI,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(list(
    summary_stats = summary_stats,
    empirical_pvalues = empirical_pvalues,
    statistical_tests = statistical_tests,
    pattern_analysis = pattern_analysis
  ))
}

#' Run Result 3 analysis for a single tree
#'
#' @param dist_obj Distance object
#' @param subset_size Subset size
#' @param n_reps Number of simulation replicates
#' @param n_random_subsets Number of random subsets
#' @param ou_alpha_values Vector of alpha values for OU process
#' @return Comprehensive Result 3 analysis results
run_result3_for_tree <- function(dist_obj, subset_size, n_reps = 100, 
                                 n_random_subsets = 100, ou_alpha_values = c(0.2, 1, 5)) {
  
  tree_name <- ifelse("tree_name" %in% names(dist_obj), 
                     dist_obj$tree_name, 
                     paste0("tree_", length(dist_obj$tip_labels), "tips"))
  
  cat(paste0("\n",strrep("=", 60), "\n"))
  cat("Result 3 analysis for", tree_name, "\n")
  cat("Subset size:", subset_size, "\n")
  cat("Simulation replicates:", n_reps, "\n")
  cat("Random subsets:", n_random_subsets, "\n")
  cat(paste0(strrep("=", 60), "\n"))
  
  # Step 1: Prepare subsets
  cat("Step 1: Preparing subsets...\n")
  subset_data <- prepare_subsets_for_trait_analysis(dist_obj, subset_size, n_random_subsets)
  
  # Step 2: Run trait analysis
  cat("\nStep 2: Running trait simulation and analysis...\n")
  trait_analysis <- run_trait_analysis(dist_obj$tree, subset_data, n_reps, ou_alpha_values)
  
  # Step 3: Analyze results
  cat("\nStep 3: Analyzing trait signals...\n")
  analysis_results <- analyze_trait_signals(trait_analysis)
  
  # Step 4: Print results
  cat("\nStep 4: Results summary\n")
  cat(paste0(strrep("-", 60), "\n"))
  
  # Print summary statistics
  if (!is.null(analysis_results$summary_stats)) {
    cat("Summary statistics (mean ± sd):\n\n")
    
    for (model in unique(analysis_results$summary_stats$Model)) {
      cat("Model:", model, "\n")
      model_stats <- analysis_results$summary_stats[analysis_results$summary_stats$Model == model, ]
      
      for (i in 1:nrow(model_stats)) {
        row <- model_stats[i, ]
        cat(sprintf("  %-10s: K = %.3f ± %.3f, λ = %.3f ± %.3f, Moran's I = %.3f ± %.3f (n=%d)\n",
                    row$Subset_Type, row$K_mean, row$K_sd, 
                    row$Lambda_mean, row$Lambda_sd,
                    row$MoransI_mean, row$MoransI_sd, row$N))
      }
      cat("\n")
    }
  }
  
  # Print pattern analysis
  if (!is.null(analysis_results$pattern_analysis) && nrow(analysis_results$pattern_analysis) > 0) {
    cat("Pattern analysis (clustered > random > dispersed):\n\n")
    
    for (i in 1:nrow(analysis_results$pattern_analysis)) {
      row <- analysis_results$pattern_analysis[i, ]
      cat("Model:", row$Model, "\n")
      cat("  K pattern holds:", row$Pattern_K, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_K), 
          "> Random:", sprintf("%.3f", row$Random_Mean_K),
          "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_K), "\n")
      cat("  Lambda pattern holds:", row$Pattern_Lambda, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_Lambda), 
          "> Random:", sprintf("%.3f", row$Random_Mean_Lambda),
          "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_Lambda), "\n")
      cat("  Moran's I pattern holds:", row$Pattern_MoransI, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_MoransI), 
          "> Random:", sprintf("%.3f", row$Random_Mean_MoransI),
          "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_MoransI), "\n")
      cat("\n")
    }
  }
  
  # Print empirical p-values
  if (!is.null(analysis_results$empirical_pvalues) && nrow(analysis_results$empirical_pvalues) > 0) {
    cat("Empirical p-values (dispersed < clustered/random):\n\n")
    
    for (model in unique(analysis_results$empirical_pvalues$Model)) {
      cat("Model:", model, "\n")
      model_pvals <- analysis_results$empirical_pvalues[analysis_results$empirical_pvalues$Model == model, ]
      
      for (metric in c("K", "Lambda", "MoransI")) {
        metric_pvals <- model_pvals[model_pvals$Metric == metric, ]
        if (nrow(metric_pvals) > 0) {
          row <- metric_pvals[1, ]
          cat(sprintf("  %-10s: dispersed (%.3f) vs clustered (%.3f): p = %.3f, vs random (%.3f): p = %.3f\n",
                      metric, row$Mean_Dispersed, row$Mean_Clustered, 
                      row$P_Dispersed_vs_Clustered, row$Mean_Random, row$P_Dispersed_vs_Random))
        }
      }
      cat("\n")
    }
  }
  
  # Print statistical test results
  if (!is.null(analysis_results$statistical_tests)) {
    cat("Statistical tests:\n\n")
    
    for (test_name in names(analysis_results$statistical_tests)) {
      if (grepl("Kruskal_Wallis", test_name)) {
        test <- analysis_results$statistical_tests[[test_name]]
        cat(test_name, ":\n")
        cat("  K: p =", test$K$p.value, "\n")
        cat("  Lambda: p =", test$Lambda$p.value, "\n")
        cat("  Moran's I: p =", test$MoransI$p.value, "\n")
        cat("\n")
      }
    }
  }
  
  cat(paste0("\n",strrep("-", 60), "\n"))
  
  # Conclusions
  cat("\nStep 5: Conclusions\n")
  cat(paste0("\n",strrep("-", 60), "\n"))
  
  if (!is.null(analysis_results$pattern_analysis) && nrow(analysis_results$pattern_analysis) > 0) {
    pattern_holds_K <- all(analysis_results$pattern_analysis$Pattern_K)
    pattern_holds_Lambda <- all(analysis_results$pattern_analysis$Pattern_Lambda)
    pattern_holds_MoransI <- all(analysis_results$pattern_analysis$Pattern_MoransI)
    
    if (pattern_holds_K && pattern_holds_Lambda && pattern_holds_MoransI) {
      cat("The expected pattern (clustered > random > dispersed) holds for ALL models and ALL metrics\n")
      cat("This strongly supports the hypothesis that phylogenetic dispersion reduces trait phylogenetic signal\n")
    } else if (pattern_holds_K && pattern_holds_Lambda) {
      cat("The expected pattern holds for K and Lambda but not for Moran's I\n")
      cat("This provides strong support for the hypothesis from traditional signal metrics\n")
    } else if (pattern_holds_K && pattern_holds_MoransI) {
      cat("The expected pattern holds for K and Moran's I but not for Lambda\n")
      cat("This provides strong support for the hypothesis from K and phylogenetic autocorrelation\n")
    } else if (pattern_holds_Lambda && pattern_holds_MoransI) {
      cat("The expected pattern holds for Lambda and Moran's I but not for K\n")
      cat("This provides strong support for the hypothesis from Lambda and phylogenetic autocorrelation\n")
    } else if (pattern_holds_K) {
      cat("The expected pattern holds for K only\n")
      cat("This provides partial support for the hypothesis\n")
    } else if (pattern_holds_Lambda) {
      cat("The expected pattern holds for Lambda only\n")
      cat("This provides partial support for the hypothesis\n")
    } else if (pattern_holds_MoransI) {
      cat("The expected pattern holds for Moran's I only\n")
      cat("This provides partial support for the hypothesis from phylogenetic autocorrelation perspective\n")
    } else {
      cat("The expected pattern does NOT hold consistently for any metric\n")
      cat("This does not support the hypothesis that phylogenetic dispersion reduces trait phylogenetic signal\n")
    }
  } else {
    cat("Insufficient data for pattern analysis\n")
  }
  
  return(list(
    subset_data = subset_data,
    trait_analysis = trait_analysis,
    analysis_results = analysis_results
  ))
}

#' Run Result 3 analysis for all trees
#'
#' @param dist_objs List of distance objects for all trees
#' @param cfg Configuration list
#' @return Comprehensive Result 3 results for all trees
run_result3_analysis <- function(dist_objs, cfg) {
  cat(paste0("\n",strrep("=", 70), "\n"))
  cat("RESULT 3: Trait-level consequences of phylogenetic dispersion\n")
  cat(paste0("\n",strrep("=", 70), "\n"))
  
  results <- list()
  
  # Process large balanced tree
  if ("large_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing large balanced tree (n =", cfg$large_n, ")\n")
    result_lb <- run_result3_for_tree(
      dist_objs$large_balanced,
      subset_size = cfg$subset_large,
      n_reps = min(100, cfg$trait_reps),  # Use smaller number for testing
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
    results$large_balanced <- result_lb
  }
  
  # Process large ladder tree
  if ("large_ladder" %in% names(dist_objs)) {
    cat("\n>>> Processing large ladder tree (n =", cfg$large_n, ")\n")
    result_ll <- run_result3_for_tree(
      dist_objs$large_ladder,
      subset_size = cfg$subset_large,
      n_reps = min(100, cfg$trait_reps),
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
    results$large_ladder <- result_ll
  }
  
  # Process small balanced tree
  if ("small_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing small balanced tree (n =", cfg$small_n, ")\n")
    result_sb <- run_result3_for_tree(
      dist_objs$small_balanced,
      subset_size = cfg$subset_small,
      n_reps = min(100, cfg$trait_reps),
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
    results$small_balanced <- result_sb
  }
  
  # Process small ladder tree
  if ("small_ladder" %in% names(dist_objs)) {
    cat("\n>>> Processing small ladder tree (n =", cfg$small_n, ")\n")
    result_sl <- run_result3_for_tree(
      dist_objs$small_ladder,
      subset_size = cfg$subset_small,
      n_reps = min(100, cfg$trait_reps),
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
    results$small_ladder <- result_sl
  }
  
  # Create overall summary
  overall_summary <- data.frame()
  
  for (tree_name in names(results)) {
    result <- results[[tree_name]]
    
    if (!is.null(result$analysis_results$pattern_analysis)) {
      pattern_df <- result$analysis_results$pattern_analysis
      
      for (i in 1:nrow(pattern_df)) {
        overall_summary <- rbind(overall_summary, data.frame(
          Tree = tree_name,
          Model = pattern_df$Model[i],
          Pattern_K = pattern_df$Pattern_K[i],
          Pattern_Lambda = pattern_df$Pattern_Lambda[i],
          Pattern_MoransI = pattern_df$Pattern_MoransI[i],
          Clustered_Mean_K = pattern_df$Clustered_Mean_K[i],
          Random_Mean_K = pattern_df$Random_Mean_K[i],
          Dispersed_Mean_K = pattern_df$Dispersed_Mean_K[i],
          Clustered_Mean_Lambda = pattern_df$Clustered_Mean_Lambda[i],
          Random_Mean_Lambda = pattern_df$Random_Mean_Lambda[i],
          Dispersed_Mean_Lambda = pattern_df$Dispersed_Mean_Lambda[i],
          Clustered_Mean_MoransI = pattern_df$Clustered_Mean_MoransI[i],
          Random_Mean_MoransI = pattern_df$Random_Mean_MoransI[i],
          Dispersed_Mean_MoransI = pattern_df$Dispersed_Mean_MoransI[i],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  cat(paste0("\n",strrep("=", 70), "\n"))
  cat("OVERALL SUMMARY\n")
  cat(paste0("\n",strrep("=", 70), "\n"))
  
  if (nrow(overall_summary) > 0) {
    print(overall_summary)
    
    # Count how many patterns hold
    n_models <- nrow(overall_summary)
    n_pattern_K <- sum(overall_summary$Pattern_K)
    n_pattern_Lambda <- sum(overall_summary$Pattern_Lambda)
    n_pattern_MoransI <- sum(overall_summary$Pattern_MoransI)
    
    cat("\nPattern summary:\n")
    cat("  K pattern holds:", n_pattern_K, "out of", n_models, "models (", 
        round(n_pattern_K/n_models*100, 1), "%)\n")
    cat("  Lambda pattern holds:", n_pattern_Lambda, "out of", n_models, "models (", 
        round(n_pattern_Lambda/n_models*100, 1), "%)\n")
    cat("  Moran's I pattern holds:", n_pattern_MoransI, "out of", n_models, "models (", 
        round(n_pattern_MoransI/n_models*100, 1), "%)\n")
    cat("  All three patterns hold:", sum(overall_summary$Pattern_K & overall_summary$Pattern_Lambda & overall_summary$Pattern_MoransI), 
        "out of", n_models, "models\n")
    cat("  K and Lambda patterns hold:", sum(overall_summary$Pattern_K & overall_summary$Pattern_Lambda), 
        "out of", n_models, "models\n")
    cat("  K and Moran's I patterns hold:", sum(overall_summary$Pattern_K & overall_summary$Pattern_MoransI), 
        "out of", n_models, "models\n")
    cat("  Lambda and Moran's I patterns hold:", sum(overall_summary$Pattern_Lambda & overall_summary$Pattern_MoransI), 
        "out of", n_models, "models\n")
  } else {
    cat("No pattern analysis results available\n")
  }
  
  return(list(
    results = results,
    overall_summary = overall_summary
  ))
}

#' Save Result 3 results to files
#'
#' @param result3_results Results from run_result3_analysis
#' @param cfg Configuration list
save_result3 <- function(result3_results, cfg) {
  cat("\nSaving Result 3 results...\n")
  
  # Create output directory if it doesn't exist
  output_dir <- cfg$tables_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save overall summary
  if (nrow(result3_results$overall_summary) > 0) {
    summary_file <- file.path(output_dir, "result3_summary.csv")
    write.csv(result3_results$overall_summary, summary_file, row.names = FALSE)
    cat("  Saved overall summary to:", summary_file, "\n")
  }
  
  # Save detailed results for each tree
  for (tree_name in names(result3_results$results)) {
    result <- result3_results$results[[tree_name]]
    
    # Save signal metrics
    if (!is.null(result$trait_analysis$signal_results)) {
      signal_file <- file.path(output_dir, paste0("result3_signal_", tree_name, ".csv"))
      write.csv(result$trait_analysis$signal_results, signal_file, row.names = FALSE)
    }
    
    # Save summary statistics
    if (!is.null(result$analysis_results$summary_stats)) {
      stats_file <- file.path(output_dir, paste0("result3_stats_", tree_name, ".csv"))
      write.csv(result$analysis_results$summary_stats, stats_file, row.names = FALSE)
    }
    
    # Save pattern analysis
    if (!is.null(result$analysis_results$pattern_analysis)) {
      pattern_file <- file.path(output_dir, paste0("result3_pattern_", tree_name, ".csv"))
      write.csv(result$analysis_results$pattern_analysis, pattern_file, row.names = FALSE)
    }
    
    # Save empirical p-values
    if (!is.null(result$analysis_results$empirical_pvalues)) {
      pvalues_file <- file.path(output_dir, paste0("result3_pvalues_", tree_name, ".csv"))
      write.csv(result$analysis_results$empirical_pvalues, pvalues_file, row.names = FALSE)
      cat("  Saved empirical p-values to:", pvalues_file, "\n")
    }
  }
  
  cat("Result 3 results saved to:", output_dir, "\n")
}

#' Test Result 3 analysis
test_result3 <- function() {
  # Load required libraries and functions
  library(ape)
  source("config/analysis_config.R")
  source("R/distance_metrics.R")
  source("R/objective_compare.R")
  source("R/subset_greedy.R")
  source("R/subset_exchange.R")
  source("R/subset_random.R")
  source("R/trait_simulation.R")
  source("R/signal_metrics.R")
  
  # Create a test configuration
  test_cfg <- list(
    large_n = 30,  # Smaller for testing
    small_n = 15,
    subset_large = 5,
    subset_small = 3,
    trait_reps = 10,
    random_subset_reps_for_trait = 20,
    ou_alpha = c(0.2, 1)
  )
  
  # Create a test tree
  cat("Creating test tree...\n")
  test_tree <- rtree(test_cfg$large_n)
  test_tree$tip.label <- paste0("sp", 1:test_cfg$large_n)
  dist_obj <- create_distance_object(test_tree)
  dist_obj$tree_name <- "test_tree"
  
  # Run Result 3 analysis
  cat("\nRunning Result 3 analysis...\n")
  result3 <- run_result3_for_tree(
    dist_obj, 
    subset_size = test_cfg$subset_large,
    n_reps = test_cfg$trait_reps,
    n_random_subsets = test_cfg$random_subset_reps_for_trait,
    ou_alpha_values = test_cfg$ou_alpha
  )
  
  # Save results
  test_cfg$tables_dir <- "test_output"
  if (!dir.exists(test_cfg$tables_dir)) {
    dir.create(test_cfg$tables_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create a simple result structure for saving
  simple_result <- list(
    results = list(test_tree = result3),
    overall_summary = if (!is.null(result3$analysis_results$pattern_analysis)) {
      result3$analysis_results$pattern_analysis
    } else {
      data.frame()
    }
  )
  
  save_result3(simple_result, test_cfg)
  
  return(result3)
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_result <- test_result3()
  cat("\nTest completed successfully.\n")
}
