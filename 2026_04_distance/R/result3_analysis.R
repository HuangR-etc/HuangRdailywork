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
#' @param metric_name Metric name to analyze
#' @return Data frame with empirical p-values for comparisons
calculate_empirical_pvalues <- function(signal_df, model_name, metric_name) {
  model_data <- signal_df[signal_df$Model == model_name, ]

  if (nrow(model_data) == 0 || !(metric_name %in% colnames(model_data))) {
    return(NULL)
  }

  tol <- 1e-15
  ess_metrics <- c("ess_1", "ess_2")

  dispersed_vals <- model_data[[metric_name]][model_data$Subset_Type == "dispersed"]
  clustered_vals <- model_data[[metric_name]][model_data$Subset_Type == "clustered"]
  random_vals <- model_data[[metric_name]][model_data$Subset_Type == "random"]

  dispersed_vals <- dispersed_vals[!is.na(dispersed_vals)]
  clustered_vals <- clustered_vals[!is.na(clustered_vals)]
  random_vals <- random_vals[!is.na(random_vals)]

  count_extreme_one_vs_many <- function(target_vals, reference_vals, direction = c("less", "greater")) {
    direction <- match.arg(direction)

    if (length(target_vals) == 0 || length(reference_vals) == 0) {
      return(NA_real_)
    }

    count_extreme <- 0
    for (target_val in target_vals) {
      if (direction == "less") {
        count_extreme <- count_extreme + sum(reference_vals >= target_val | abs(reference_vals - target_val) < tol)
      } else {
        count_extreme <- count_extreme + sum(reference_vals <= target_val | abs(reference_vals - target_val) < tol)
      }
    }

    count_extreme <- count_extreme / length(target_vals)
    (count_extreme + 1) / (length(reference_vals) + 1)
  }

  count_extreme_pairwise <- function(left_vals, right_vals, direction = c("less", "greater")) {
    direction <- match.arg(direction)

    if (length(left_vals) == 0 || length(right_vals) == 0) {
      return(NA_real_)
    }

    total_comparisons <- length(left_vals) * length(right_vals)
    count_extreme <- 0

    for (left_val in left_vals) {
      if (direction == "less") {
        count_extreme <- count_extreme + sum(right_vals >= left_val | abs(right_vals - left_val) < tol)
      } else {
        count_extreme <- count_extreme + sum(right_vals <= left_val | abs(right_vals - left_val) < tol)
      }
    }

    (count_extreme + 1) / (total_comparisons + 1)
  }

  if (metric_name == "MoransI") {
    dispersed_abs <- abs(dispersed_vals)
    clustered_abs <- abs(clustered_vals)
    random_abs <- abs(random_vals)

    p_dispersed_vs_clustered <- count_extreme_one_vs_many(dispersed_abs, clustered_abs, direction = "less")
    p_dispersed_vs_random <- count_extreme_one_vs_many(dispersed_abs, random_abs, direction = "less")
    p_random_vs_clustered <- count_extreme_pairwise(random_abs, clustered_abs, direction = "less")
  } else if (metric_name %in% ess_metrics) {
    p_dispersed_vs_clustered <- count_extreme_one_vs_many(dispersed_vals, clustered_vals, direction = "greater")
    p_dispersed_vs_random <- count_extreme_one_vs_many(dispersed_vals, random_vals, direction = "greater")
    p_random_vs_clustered <- count_extreme_pairwise(random_vals, clustered_vals, direction = "greater")
  } else {
    p_dispersed_vs_clustered <- count_extreme_one_vs_many(dispersed_vals, clustered_vals, direction = "less")
    p_dispersed_vs_random <- count_extreme_one_vs_many(dispersed_vals, random_vals, direction = "less")
    p_random_vs_clustered <- count_extreme_pairwise(random_vals, clustered_vals, direction = "less")
  }

  mean_dispersed <- mean(dispersed_vals, na.rm = TRUE)
  mean_clustered <- mean(clustered_vals, na.rm = TRUE)
  mean_random <- mean(random_vals, na.rm = TRUE)

  data.frame(
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
  )
}

#' Analyze trait signal results
#'
#' @param trait_analysis_results Results from run_trait_analysis
#' @return Statistical analysis of trait signals
analyze_trait_signals <- function(trait_analysis_results) {
  cat("Analyzing trait signal results...\n")

  signal_df <- trait_analysis_results$signal_results
  metric_names <- c(
    "K", "Lambda", "MoransI", "mpnns",
    "mean_offdiag_cor", "max_offdiag_cor", "ess_1", "ess_2"
  )
  ess_metrics <- c("ess_1", "ess_2")

  signal_df <- signal_df[complete.cases(signal_df[, metric_names]), ]

  if (nrow(signal_df) == 0) {
    warning("No complete cases found in signal results")
    return(NULL)
  }

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
          mpnns_mean = mean(subset_data$mpnns, na.rm = TRUE),
          mpnns_sd = sd(subset_data$mpnns, na.rm = TRUE),
          mpnns_min = min(subset_data$mpnns, na.rm = TRUE),
          mpnns_max = max(subset_data$mpnns, na.rm = TRUE),
          mean_offdiag_cor_mean = mean(subset_data$mean_offdiag_cor, na.rm = TRUE),
          mean_offdiag_cor_sd = sd(subset_data$mean_offdiag_cor, na.rm = TRUE),
          mean_offdiag_cor_min = min(subset_data$mean_offdiag_cor, na.rm = TRUE),
          mean_offdiag_cor_max = max(subset_data$mean_offdiag_cor, na.rm = TRUE),
          max_offdiag_cor_mean = mean(subset_data$max_offdiag_cor, na.rm = TRUE),
          max_offdiag_cor_sd = sd(subset_data$max_offdiag_cor, na.rm = TRUE),
          max_offdiag_cor_min = min(subset_data$max_offdiag_cor, na.rm = TRUE),
          max_offdiag_cor_max = max(subset_data$max_offdiag_cor, na.rm = TRUE),
          ess_1_mean = mean(subset_data$ess_1, na.rm = TRUE),
          ess_1_sd = sd(subset_data$ess_1, na.rm = TRUE),
          ess_1_min = min(subset_data$ess_1, na.rm = TRUE),
          ess_1_max = max(subset_data$ess_1, na.rm = TRUE),
          ess_2_mean = mean(subset_data$ess_2, na.rm = TRUE),
          ess_2_sd = sd(subset_data$ess_2, na.rm = TRUE),
          ess_2_min = min(subset_data$ess_2, na.rm = TRUE),
          ess_2_max = max(subset_data$ess_2, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  empirical_pvalues <- data.frame()
  for (model in models) {
    for (metric in metric_names) {
      pval_df <- calculate_empirical_pvalues(signal_df, model, metric)
      if (!is.null(pval_df)) {
        empirical_pvalues <- rbind(empirical_pvalues, pval_df)
      }
    }
  }

  statistical_tests <- list()
  for (model in models) {
    model_data <- signal_df[signal_df$Model == model, ]

    if (nrow(model_data) > 0) {
      if (length(unique(model_data$Subset_Type)) >= 2) {
        statistical_tests[[paste0(model, "_Kruskal_Wallis")]] <- list(
          K = kruskal.test(K ~ Subset_Type, data = model_data),
          Lambda = kruskal.test(Lambda ~ Subset_Type, data = model_data),
          MoransI = kruskal.test(MoransI ~ Subset_Type, data = model_data),
          mpnns = kruskal.test(mpnns ~ Subset_Type, data = model_data),
          mean_offdiag_cor = kruskal.test(mean_offdiag_cor ~ Subset_Type, data = model_data),
          max_offdiag_cor = kruskal.test(max_offdiag_cor ~ Subset_Type, data = model_data),
          ess_1 = kruskal.test(ess_1 ~ Subset_Type, data = model_data),
          ess_2 = kruskal.test(ess_2 ~ Subset_Type, data = model_data)
        )
      }

      pairwise_results <- list()
      subset_pairs <- list(
        c("clustered", "random"),
        c("random", "dispersed"),
        c("clustered", "dispersed")
      )

      for (pair in subset_pairs) {
        type1 <- pair[1]
        type2 <- pair[2]
        pair_key <- paste0(type1, "_vs_", type2)
        pairwise_results[[pair_key]] <- list()

        for (metric in metric_names) {
          data1 <- model_data[[metric]][model_data$Subset_Type == type1]
          data2 <- model_data[[metric]][model_data$Subset_Type == type2]
          data1 <- data1[!is.na(data1)]
          data2 <- data2[!is.na(data2)]

          if (length(data1) > 1 && length(data2) > 1) {
            alt <- if (metric %in% ess_metrics) "less" else "greater"
            pairwise_results[[pair_key]][[metric]] <- wilcox.test(data1, data2, alternative = alt)
          } else {
            pairwise_results[[pair_key]][[metric]] <- NULL
          }
        }
      }

      statistical_tests[[paste0(model, "_Pairwise")]] <- pairwise_results
    }
  }

  pattern_analysis <- data.frame()
  tolerance <- 1e-10

  for (model in models) {
    model_summary <- summary_stats[summary_stats$Model == model, ]

    if (nrow(model_summary) >= 3) {
      clustered_mean_K <- model_summary$K_mean[model_summary$Subset_Type == "clustered"]
      random_mean_K <- model_summary$K_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_K <- model_summary$K_mean[model_summary$Subset_Type == "dispersed"]

      clustered_mean_Lambda <- model_summary$Lambda_mean[model_summary$Subset_Type == "clustered"]
      random_mean_Lambda <- model_summary$Lambda_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_Lambda <- model_summary$Lambda_mean[model_summary$Subset_Type == "dispersed"]

      clustered_mean_MoransI <- model_summary$MoransI_mean[model_summary$Subset_Type == "clustered"]
      random_mean_MoransI <- model_summary$MoransI_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_MoransI <- model_summary$MoransI_mean[model_summary$Subset_Type == "dispersed"]

      clustered_mean_mpnns <- model_summary$mpnns_mean[model_summary$Subset_Type == "clustered"]
      random_mean_mpnns <- model_summary$mpnns_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_mpnns <- model_summary$mpnns_mean[model_summary$Subset_Type == "dispersed"]

      clustered_mean_mean_offdiag_cor <- model_summary$mean_offdiag_cor_mean[model_summary$Subset_Type == "clustered"]
      random_mean_mean_offdiag_cor <- model_summary$mean_offdiag_cor_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_mean_offdiag_cor <- model_summary$mean_offdiag_cor_mean[model_summary$Subset_Type == "dispersed"]

      clustered_mean_max_offdiag_cor <- model_summary$max_offdiag_cor_mean[model_summary$Subset_Type == "clustered"]
      random_mean_max_offdiag_cor <- model_summary$max_offdiag_cor_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_max_offdiag_cor <- model_summary$max_offdiag_cor_mean[model_summary$Subset_Type == "dispersed"]

      clustered_mean_ess_1 <- model_summary$ess_1_mean[model_summary$Subset_Type == "clustered"]
      random_mean_ess_1 <- model_summary$ess_1_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_ess_1 <- model_summary$ess_1_mean[model_summary$Subset_Type == "dispersed"]

      clustered_mean_ess_2 <- model_summary$ess_2_mean[model_summary$Subset_Type == "clustered"]
      random_mean_ess_2 <- model_summary$ess_2_mean[model_summary$Subset_Type == "random"]
      dispersed_mean_ess_2 <- model_summary$ess_2_mean[model_summary$Subset_Type == "dispersed"]

      pattern_analysis <- rbind(pattern_analysis, data.frame(
        Model = model,
        Pattern_K = (clustered_mean_K > random_mean_K - tolerance) && (random_mean_K > dispersed_mean_K - tolerance),
        Pattern_Lambda = (clustered_mean_Lambda > random_mean_Lambda - tolerance) && (random_mean_Lambda > dispersed_mean_Lambda - tolerance),
        Pattern_MoransI = (clustered_mean_MoransI > random_mean_MoransI - tolerance) && (random_mean_MoransI > dispersed_mean_MoransI - tolerance),
        Pattern_mpnns = (clustered_mean_mpnns > random_mean_mpnns - tolerance) && (random_mean_mpnns > dispersed_mean_mpnns - tolerance),
        Pattern_mean_offdiag_cor = (clustered_mean_mean_offdiag_cor > random_mean_mean_offdiag_cor - tolerance) && (random_mean_mean_offdiag_cor > dispersed_mean_mean_offdiag_cor - tolerance),
        Pattern_max_offdiag_cor = (clustered_mean_max_offdiag_cor > random_mean_max_offdiag_cor - tolerance) && (random_mean_max_offdiag_cor > dispersed_mean_max_offdiag_cor - tolerance),
        Pattern_ess_1 = (dispersed_mean_ess_1 > random_mean_ess_1 - tolerance) && (random_mean_ess_1 > clustered_mean_ess_1 - tolerance),
        Pattern_ess_2 = (dispersed_mean_ess_2 > random_mean_ess_2 - tolerance) && (random_mean_ess_2 > clustered_mean_ess_2 - tolerance),
        Clustered_Mean_K = clustered_mean_K,
        Random_Mean_K = random_mean_K,
        Dispersed_Mean_K = dispersed_mean_K,
        Clustered_Mean_Lambda = clustered_mean_Lambda,
        Random_Mean_Lambda = random_mean_Lambda,
        Dispersed_Mean_Lambda = dispersed_mean_Lambda,
        Clustered_Mean_MoransI = clustered_mean_MoransI,
        Random_Mean_MoransI = random_mean_MoransI,
        Dispersed_Mean_MoransI = dispersed_mean_MoransI,
        Clustered_Mean_mpnns = clustered_mean_mpnns,
        Random_Mean_mpnns = random_mean_mpnns,
        Dispersed_Mean_mpnns = dispersed_mean_mpnns,
        Clustered_Mean_mean_offdiag_cor = clustered_mean_mean_offdiag_cor,
        Random_Mean_mean_offdiag_cor = random_mean_mean_offdiag_cor,
        Dispersed_Mean_mean_offdiag_cor = dispersed_mean_mean_offdiag_cor,
        Clustered_Mean_max_offdiag_cor = clustered_mean_max_offdiag_cor,
        Random_Mean_max_offdiag_cor = random_mean_max_offdiag_cor,
        Dispersed_Mean_max_offdiag_cor = dispersed_mean_max_offdiag_cor,
        Clustered_Mean_ess_1 = clustered_mean_ess_1,
        Random_Mean_ess_1 = random_mean_ess_1,
        Dispersed_Mean_ess_1 = dispersed_mean_ess_1,
        Clustered_Mean_ess_2 = clustered_mean_ess_2,
        Random_Mean_ess_2 = random_mean_ess_2,
        Dispersed_Mean_ess_2 = dispersed_mean_ess_2,
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

  cat(paste0("\n", strrep("=", 60), "\n"))
  cat("Result 3 analysis for", tree_name, "\n")
  cat("Subset size:", subset_size, "\n")
  cat("Simulation replicates:", n_reps, "\n")
  cat("Random subsets:", n_random_subsets, "\n")
  cat(paste0(strrep("=", 60), "\n"))

  cat("Step 1: Preparing subsets...\n")
  subset_data <- prepare_subsets_for_trait_analysis(dist_obj, subset_size, n_random_subsets)

  cat("\nStep 2: Running trait simulation and analysis...\n")
  trait_analysis <- run_trait_analysis(dist_obj$tree, subset_data, n_reps, ou_alpha_values)

  cat("\nStep 3: Analyzing trait signals...\n")
  analysis_results <- analyze_trait_signals(trait_analysis)

  cat("\nStep 4: Results summary\n")
  cat(paste0(strrep("-", 60), "\n"))

  if (!is.null(analysis_results$summary_stats)) {
    cat("Summary statistics (mean ± sd):\n\n")
    for (model in unique(analysis_results$summary_stats$Model)) {
      cat("Model:", model, "\n")
      model_stats <- analysis_results$summary_stats[analysis_results$summary_stats$Model == model, ]

      for (i in 1:nrow(model_stats)) {
        row <- model_stats[i, ]
        cat(sprintf(
          "  %-10s: K = %.3f ± %.3f, λ = %.3f ± %.3f, Moran's I = %.3f ± %.3f, mpnns = %.3f ± %.3f (n=%d)\n",
          row$Subset_Type, row$K_mean, row$K_sd,
          row$Lambda_mean, row$Lambda_sd,
          row$MoransI_mean, row$MoransI_sd,
          row$mpnns_mean, row$mpnns_sd, row$N
        ))
        cat(sprintf(
          "              mean_offdiag_cor = %.3f ± %.3f, max_offdiag_cor = %.3f ± %.3f, ess_1 = %.3f ± %.3f, ess_2 = %.3f ± %.3f\n",
          row$mean_offdiag_cor_mean, row$mean_offdiag_cor_sd,
          row$max_offdiag_cor_mean, row$max_offdiag_cor_sd,
          row$ess_1_mean, row$ess_1_sd,
          row$ess_2_mean, row$ess_2_sd
        ))
      }
      cat("\n")
    }
  }

  if (!is.null(analysis_results$pattern_analysis) && nrow(analysis_results$pattern_analysis) > 0) {
    cat("Pattern analysis:\n")
    cat("  Expected directions: clustered > random > dispersed for K, Lambda, Moran's I, mpnns, mean_offdiag_cor, max_offdiag_cor\n")
    cat("                       dispersed > random > clustered for ess_1 and ess_2\n\n")

    for (i in 1:nrow(analysis_results$pattern_analysis)) {
      row <- analysis_results$pattern_analysis[i, ]
      cat("Model:", row$Model, "\n")
      cat("  K pattern holds:", row$Pattern_K, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_K), "> Random:", sprintf("%.3f", row$Random_Mean_K), "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_K), "\n")
      cat("  Lambda pattern holds:", row$Pattern_Lambda, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_Lambda), "> Random:", sprintf("%.3f", row$Random_Mean_Lambda), "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_Lambda), "\n")
      cat("  Moran's I pattern holds:", row$Pattern_MoransI, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_MoransI), "> Random:", sprintf("%.3f", row$Random_Mean_MoransI), "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_MoransI), "\n")
      cat("  mpnns pattern holds:", row$Pattern_mpnns, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_mpnns), "> Random:", sprintf("%.3f", row$Random_Mean_mpnns), "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_mpnns), "\n")
      cat("  mean_offdiag_cor pattern holds:", row$Pattern_mean_offdiag_cor, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_mean_offdiag_cor), "> Random:", sprintf("%.3f", row$Random_Mean_mean_offdiag_cor), "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_mean_offdiag_cor), "\n")
      cat("  max_offdiag_cor pattern holds:", row$Pattern_max_offdiag_cor, "\n")
      cat("    Clustered:", sprintf("%.3f", row$Clustered_Mean_max_offdiag_cor), "> Random:", sprintf("%.3f", row$Random_Mean_max_offdiag_cor), "> Dispersed:", sprintf("%.3f", row$Dispersed_Mean_max_offdiag_cor), "\n")
      cat("  ess_1 pattern holds:", row$Pattern_ess_1, "\n")
      cat("    Dispersed:", sprintf("%.3f", row$Dispersed_Mean_ess_1), "> Random:", sprintf("%.3f", row$Random_Mean_ess_1), "> Clustered:", sprintf("%.3f", row$Clustered_Mean_ess_1), "\n")
      cat("  ess_2 pattern holds:", row$Pattern_ess_2, "\n")
      cat("    Dispersed:", sprintf("%.3f", row$Dispersed_Mean_ess_2), "> Random:", sprintf("%.3f", row$Random_Mean_ess_2), "> Clustered:", sprintf("%.3f", row$Clustered_Mean_ess_2), "\n\n")
    }
  }

  if (!is.null(analysis_results$empirical_pvalues) && nrow(analysis_results$empirical_pvalues) > 0) {
    cat("Empirical p-values (expected one-sided direction by metric):\n\n")
    metric_order <- c("K", "Lambda", "MoransI", "mpnns", "mean_offdiag_cor", "max_offdiag_cor", "ess_1", "ess_2")

    for (model in unique(analysis_results$empirical_pvalues$Model)) {
      cat("Model:", model, "\n")
      model_pvals <- analysis_results$empirical_pvalues[analysis_results$empirical_pvalues$Model == model, ]

      for (metric in metric_order) {
        metric_pvals <- model_pvals[model_pvals$Metric == metric, ]
        if (nrow(metric_pvals) > 0) {
          row <- metric_pvals[1, ]
          cat(sprintf(
            "  %-16s: dispersed (%.3f) vs clustered (%.3f): p = %.3f, vs random (%.3f): p = %.3f\n",
            metric, row$Mean_Dispersed, row$Mean_Clustered,
            row$P_Dispersed_vs_Clustered, row$Mean_Random, row$P_Dispersed_vs_Random
          ))
        }
      }
      cat("\n")
    }
  }

  if (!is.null(analysis_results$statistical_tests)) {
    cat("Statistical tests:\n\n")
    metric_order <- c("K", "Lambda", "MoransI", "mpnns", "mean_offdiag_cor", "max_offdiag_cor", "ess_1", "ess_2")

    for (test_name in names(analysis_results$statistical_tests)) {
      if (grepl("Kruskal_Wallis", test_name)) {
        test <- analysis_results$statistical_tests[[test_name]]
        cat(test_name, ":\n")
        for (metric in metric_order) {
          if (!is.null(test[[metric]])) {
            cat(" ", sprintf("%-16s p = %.6f", metric, test[[metric]]$p.value), "\n", sep = "")
          }
        }
        cat("\n")
      }
    }
  }

  cat(paste0("\n", strrep("-", 60), "\n"))
  cat("\nStep 5: Conclusions\n")
  cat(paste0("\n", strrep("-", 60), "\n"))

  if (!is.null(analysis_results$pattern_analysis) && nrow(analysis_results$pattern_analysis) > 0) {
    pattern_cols <- c("Pattern_K", "Pattern_Lambda", "Pattern_MoransI", "Pattern_mpnns", "Pattern_mean_offdiag_cor", "Pattern_max_offdiag_cor", "Pattern_ess_1", "Pattern_ess_2")
    pattern_holds <- vapply(pattern_cols, function(col) all(analysis_results$pattern_analysis[[col]]), logical(1))
    n_metrics_hold <- sum(pattern_holds)

    cat("Metrics with the expected ordering across all analyzed models:", n_metrics_hold, "out of", length(pattern_cols), "\n")
    for (col in names(pattern_holds)) {
      cat(" ", sprintf("%-28s %s", col, pattern_holds[[col]]), "\n", sep = "")
    }

    if (all(pattern_holds[c("Pattern_mean_offdiag_cor", "Pattern_max_offdiag_cor", "Pattern_ess_1", "Pattern_ess_2")])) {
      cat("\nThe new tree-induced dependence diagnostics consistently support the interpretation that dispersed subsets carry weaker tree-induced dependence.\n")
    } else {
      cat("\nThe new tree-induced dependence diagnostics do not all show the expected ordering, so support is mixed.\n")
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
  cat(paste0("\n", strrep("=", 70), "\n"))
  cat("RESULT 3: Trait-level consequences of phylogenetic dispersion\n")
  cat(paste0("\n", strrep("=", 70), "\n"))

  results <- list()

  if ("large_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing large balanced tree (n =", cfg$large_n, ")\n")
    results$large_balanced <- run_result3_for_tree(
      dist_objs$large_balanced,
      subset_size = cfg$subset_large,
      n_reps = cfg$trait_reps,
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
  }

  if ("large_ladder" %in% names(dist_objs)) {
    cat("\n>>> Processing large ladder tree (n =", cfg$large_n, ")\n")
    results$large_ladder <- run_result3_for_tree(
      dist_objs$large_ladder,
      subset_size = cfg$subset_large,
      n_reps = cfg$trait_reps,
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
  }

  if ("small_balanced" %in% names(dist_objs)) {
    cat("\n>>> Processing small balanced tree (n =", cfg$small_n, ")\n")
    results$small_balanced <- run_result3_for_tree(
      dist_objs$small_balanced,
      subset_size = cfg$subset_small,
      n_reps = cfg$trait_reps,
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
  }

  if ("small_ladder" %in% names(dist_objs)) {
    cat("\n>>> Processing small ladder tree (n =", cfg$small_n, ")\n")
    results$small_ladder <- run_result3_for_tree(
      dist_objs$small_ladder,
      subset_size = cfg$subset_small,
      n_reps = cfg$trait_reps,
      n_random_subsets = cfg$random_subset_reps_for_trait,
      ou_alpha_values = cfg$ou_alpha
    )
  }

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
          Pattern_mpnns = pattern_df$Pattern_mpnns[i],
          Pattern_mean_offdiag_cor = pattern_df$Pattern_mean_offdiag_cor[i],
          Pattern_max_offdiag_cor = pattern_df$Pattern_max_offdiag_cor[i],
          Pattern_ess_1 = pattern_df$Pattern_ess_1[i],
          Pattern_ess_2 = pattern_df$Pattern_ess_2[i],
          Clustered_Mean_K = pattern_df$Clustered_Mean_K[i],
          Random_Mean_K = pattern_df$Random_Mean_K[i],
          Dispersed_Mean_K = pattern_df$Dispersed_Mean_K[i],
          Clustered_Mean_Lambda = pattern_df$Clustered_Mean_Lambda[i],
          Random_Mean_Lambda = pattern_df$Random_Mean_Lambda[i],
          Dispersed_Mean_Lambda = pattern_df$Dispersed_Mean_Lambda[i],
          Clustered_Mean_MoransI = pattern_df$Clustered_Mean_MoransI[i],
          Random_Mean_MoransI = pattern_df$Random_Mean_MoransI[i],
          Dispersed_Mean_MoransI = pattern_df$Dispersed_Mean_MoransI[i],
          Clustered_Mean_mpnns = pattern_df$Clustered_Mean_mpnns[i],
          Random_Mean_mpnns = pattern_df$Random_Mean_mpnns[i],
          Dispersed_Mean_mpnns = pattern_df$Dispersed_Mean_mpnns[i],
          Clustered_Mean_mean_offdiag_cor = pattern_df$Clustered_Mean_mean_offdiag_cor[i],
          Random_Mean_mean_offdiag_cor = pattern_df$Random_Mean_mean_offdiag_cor[i],
          Dispersed_Mean_mean_offdiag_cor = pattern_df$Dispersed_Mean_mean_offdiag_cor[i],
          Clustered_Mean_max_offdiag_cor = pattern_df$Clustered_Mean_max_offdiag_cor[i],
          Random_Mean_max_offdiag_cor = pattern_df$Random_Mean_max_offdiag_cor[i],
          Dispersed_Mean_max_offdiag_cor = pattern_df$Dispersed_Mean_max_offdiag_cor[i],
          Clustered_Mean_ess_1 = pattern_df$Clustered_Mean_ess_1[i],
          Random_Mean_ess_1 = pattern_df$Random_Mean_ess_1[i],
          Dispersed_Mean_ess_1 = pattern_df$Dispersed_Mean_ess_1[i],
          Clustered_Mean_ess_2 = pattern_df$Clustered_Mean_ess_2[i],
          Random_Mean_ess_2 = pattern_df$Random_Mean_ess_2[i],
          Dispersed_Mean_ess_2 = pattern_df$Dispersed_Mean_ess_2[i],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  cat(paste0("\n", strrep("=", 70), "\n"))
  cat("OVERALL SUMMARY\n")
  cat(paste0("\n", strrep("=", 70), "\n"))

  if (nrow(overall_summary) > 0) {
    print(overall_summary)
    pattern_cols <- c("Pattern_K", "Pattern_Lambda", "Pattern_MoransI", "Pattern_mpnns", "Pattern_mean_offdiag_cor", "Pattern_max_offdiag_cor", "Pattern_ess_1", "Pattern_ess_2")
    n_models <- nrow(overall_summary)

    cat("\nPattern summary:\n")
    for (col in pattern_cols) {
      n_hold <- sum(overall_summary[[col]])
      cat(" ", sprintf("%-28s %d out of %d models (%.1f%%)", col, n_hold, n_models, 100 * n_hold / n_models), "\n", sep = "")
    }
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
  
  # Create a test configuration (using new design: 1 trait rep, 20 random subsets)
  test_cfg <- list(
    large_n = 30,  # Smaller for testing
    small_n = 15,
    subset_large = 5,
    subset_small = 3,
    trait_reps = 1,  # NEW: Only 1 trait simulation per tree
    random_subset_reps_for_trait = 20,  # 20 random subsets for testing
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
