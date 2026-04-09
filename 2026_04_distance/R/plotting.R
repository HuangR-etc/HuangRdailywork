# Plotting functions for phylogenetic dispersed subset analysis
# This module provides visualization functions for all results

library(ggplot2)
library(gridExtra)
library(viridis)

#' Create theme for publication-quality plots
#'
#' @return A ggplot2 theme
theme_phylogenetic <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold", size = 11)
    )
}

#' Plot null distribution with observed value
#'
#' @param null_metrics Data frame with null distribution metrics
#' @param observed_metrics List with observed metrics
#' @param metric_name Name of metric to plot ("MinPD", "MeanPD", or "MeanNND")
#' @param title Plot title
#' @param x_label X-axis label
#' @return A ggplot object
plot_null_distribution <- function(null_metrics, observed_metrics, metric_name, 
                                   title = NULL, x_label = NULL) {
  
  if (is.null(title)) {
    title <- paste("Null distribution of", metric_name)
  }
  
  if (is.null(x_label)) {
    x_label <- metric_name
  }
  
  observed_value <- observed_metrics[[metric_name]]
  
  p <- ggplot(null_metrics, aes(x = .data[[metric_name]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, 
                   fill = "steelblue", alpha = 0.7, color = "white") +
    geom_density(linewidth = 1, color = "darkblue") +
    geom_vline(xintercept = observed_value, color = "red", 
               linewidth = 1.5, linetype = "dashed") +
    annotate("text", x = observed_value, y = Inf, 
             label = paste("Observed:", round(observed_value, 3)),
             vjust = 2, hjust = 1.1, color = "red", size = 4) +
    labs(title = title, x = x_label, y = "Density") +
    theme_phylogenetic()
  
  return(p)
}

#' Plot comparison of multiple algorithms
#'
#' @param algorithm_metrics Data frame with metrics for multiple algorithms
#' @param metrics_to_plot Vector of metrics to plot (default: all three)
#' @param title Plot title
#' @return A list of ggplot objects
plot_algorithm_comparison <- function(algorithm_metrics, 
                                      metrics_to_plot = c("MinPD", "MeanPD", "MeanNND"),
                                      title = "Algorithm Comparison") {
  
  plots <- list()
  
  for (metric in metrics_to_plot) {
    p <- ggplot(algorithm_metrics, aes(x = Algorithm, y = .data[[metric]], fill = Algorithm)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_text(aes(label = round(.data[[metric]], 3)), vjust = -0.5, size = 3.5) +
      scale_fill_viridis(discrete = TRUE, option = "D") +
      labs(title = paste(metric, "by Algorithm"),
           x = "Algorithm", y = metric) +
      theme_phylogenetic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots[[metric]] <- p
  }
  
  # Create a combined plot
  if (length(plots) == 3) {
    combined <- grid.arrange(
      plots[[1]], plots[[2]], plots[[3]],
      ncol = 3,
      top = title
    )
    plots$combined <- combined
  }
  
  return(plots)
}

#' Plot trait signal comparison across subset types
#'
#' @param signal_results Data frame with signal metrics
#' @param metric_name Signal metric to plot ("K" or "Lambda")
#' @param title Plot title
#' @return A ggplot object
plot_trait_signal <- function(signal_results, metric_name = "K", 
                              title = "Phylogenetic Signal by Subset Type") {
  
  # Check if metric exists
  if (!metric_name %in% names(signal_results)) {
    stop("Metric", metric_name, "not found in signal results")
  }
  
  p <- ggplot(signal_results, aes(x = Subset_Type, y = .data[[metric_name]], fill = Subset_Type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    facet_wrap(~ Model, scales = "free_y") +
    scale_fill_manual(values = c("dispersed" = "#2E8B57", 
                                 "clustered" = "#CD5C5C", 
                                 "random" = "#4682B4")) +
    labs(title = title,
         x = "Subset Type", 
         y = ifelse(metric_name == "K", "Blomberg's K", "Pagel's λ")) +
    theme_phylogenetic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Plot heuristic vs exact comparison
#'
#' @param comparison_summary Data frame with comparison summary
#' @param title Plot title
#' @return A ggplot object
plot_heuristic_vs_exact <- function(comparison_summary, title = "Heuristic vs Exact Optimum") {
  
  # Reshape data for plotting
  plot_data <- data.frame(
    Metric = rep(comparison_summary$Metric, 2),
    Value = c(comparison_summary$Heuristic, comparison_summary$Exact),
    Type = rep(c("Heuristic", "Exact"), each = nrow(comparison_summary))
  )
  
  p <- ggplot(plot_data, aes(x = Metric, y = Value, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    geom_text(aes(label = round(Value, 3)), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Heuristic" = "#FF6B6B", "Exact" = "#4ECDC4")) +
    labs(title = title,
         x = "Metric", y = "Value", fill = "Algorithm") +
    theme_phylogenetic()
  
  return(p)
}

#' Plot gap between heuristic and exact optimum
#'
#' @param comparison_summary Data frame with comparison summary
#' @param title Plot title
#' @return A ggplot object
plot_gap_analysis <- function(comparison_summary, title = "Gap Analysis: Heuristic vs Exact") {
  
  p <- ggplot(comparison_summary, aes(x = Metric, y = Gap)) +
    geom_bar(stat = "identity", aes(fill = Gap > 0), alpha = 0.8) +
    geom_text(aes(label = sprintf("%+.3f\n(%+.1f%%)", Gap, Rel_Gap_Pct)), 
              vjust = ifelse(comparison_summary$Gap > 0, -0.5, 1.5), 
              size = 3.5) +
    scale_fill_manual(values = c("TRUE" = "#4ECDC4", "FALSE" = "#FF6B6B"),
                      labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                      name = "Gap Direction") +
    labs(title = title,
         x = "Metric", y = "Gap (Exact - Heuristic)") +
    theme_phylogenetic()
  
  return(p)
}

#' Create multi-panel figure for Result 1
#'
#' @param result1_results Results from Result 1 analysis
#' @param tree_name Name of the tree
#' @param output_dir Output directory for saving plots
#' @return List of plot objects
create_result1_figures <- function(result1_results, tree_name, output_dir) {
  cat("Creating Result 1 figures for", tree_name, "...\n")
  
  result <- result1_results$results[[tree_name]]
  
  plots <- list()
  
  # Plot null distributions for all three metrics
  metrics <- c("MinPD", "MeanPD", "MeanNND")
  
  for (metric in metrics) {
    p <- plot_null_distribution(
      null_metrics = result$comparison$null_metrics,
      observed_metrics = result$summary$observed_metrics,
      metric_name = metric,
      title = paste("Null Distribution of", metric, "-", tree_name)
    )
    
    plots[[paste0("null_", metric)]] <- p
    
    # Save individual plot
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, paste0("result1_null_", metric, "_", tree_name, ".pdf"))
      ggsave(plot_file, p, width = 8, height = 6)
    }
  }
  
  # Create combined null distribution plot
  if (length(plots) == 3) {
    combined_null <- grid.arrange(
      plots[[1]], plots[[2]], plots[[3]],
      ncol = 3,
      top = paste("Null Distributions -", tree_name)
    )
    plots$combined_null <- combined_null
    
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, paste0("result1_combined_null_", tree_name, ".pdf"))
      ggsave(plot_file, combined_null, width = 16, height = 6)
    }
  }
  
  # Plot z-scores
  z_scores <- result$summary$z_scores
  z_df <- data.frame(
    Metric = names(z_scores),
    Z_Score = unlist(z_scores)
  )
  
  p_z <- ggplot(z_df, aes(x = Metric, y = Z_Score, fill = Metric)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "red") +
    geom_text(aes(label = round(Z_Score, 2)), vjust = -0.5, size = 4) +
    scale_fill_viridis(discrete = TRUE, option = "C") +
    labs(title = paste("Z-Scores -", tree_name),
         x = "Metric", y = "Z-Score") +
    theme_phylogenetic()
  
  plots$z_scores <- p_z
  
  if (!is.null(output_dir)) {
    plot_file <- file.path(output_dir, paste0("result1_z_scores_", tree_name, ".pdf"))
    ggsave(plot_file, p_z, width = 8, height = 6)
  }
  
  return(plots)
}

#' Create multi-panel figure for Result 2
#'
#' @param result2_results Results from Result 2 analysis
#' @param tree_name Name of the tree
#' @param output_dir Output directory for saving plots
#' @return List of plot objects
create_result2_figures <- function(result2_results, tree_name, output_dir) {
  cat("Creating Result 2 figures for", tree_name, "...\n")
  
  result <- result2_results$results[[tree_name]]
  
  plots <- list()
  
  # Plot algorithm comparison
  algorithm_plots <- plot_algorithm_comparison(
    algorithm_metrics = result$comparison$metrics,
    title = paste("Algorithm Comparison -", tree_name)
  )
  
  plots <- c(plots, algorithm_plots)
  
  # Save individual algorithm plots
  if (!is.null(output_dir)) {
    for (metric_name in names(algorithm_plots)) {
      if (metric_name != "combined") {
        plot_file <- file.path(output_dir, paste0("result2_", metric_name, "_", tree_name, ".pdf"))
        ggsave(plot_file, algorithm_plots[[metric_name]], width = 8, height = 6)
      }
    }
    
    # Save combined plot
    if ("combined" %in% names(algorithm_plots)) {
      plot_file <- file.path(output_dir, paste0("result2_combined_", tree_name, ".pdf"))
      ggsave(plot_file, algorithm_plots$combined, width = 14, height = 6)
    }
  }
  
  # Plot improvement from greedy to exchange
  if (!is.null(result$algorithm_results$A_main$improvement)) {
    improvement <- result$algorithm_results$A_main$improvement
    imp_df <- data.frame(
      Metric = names(improvement),
      Improvement = unlist(improvement)
    )
    
    p_imp <- ggplot(imp_df, aes(x = Metric, y = Improvement, fill = Improvement > 0)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_text(aes(label = round(Improvement, 3)), vjust = -0.5, size = 4) +
      scale_fill_manual(values = c("TRUE" = "#4ECDC4", "FALSE" = "#FF6B6B"),
                        name = "Positive Improvement") +
      labs(title = paste("Improvement from Exchange Refinement -", tree_name),
           x = "Metric", y = "Improvement") +
      theme_phylogenetic()
    
    plots$improvement <- p_imp
    
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, paste0("result2_improvement_", tree_name, ".pdf"))
      ggsave(plot_file, p_imp, width = 8, height = 6)
    }
  }
  
  return(plots)
}

#' Create multi-panel figure for Result 3
#'
#' @param result3_results Results from Result 3 analysis
#' @param tree_name Name of the tree
#' @param output_dir Output directory for saving plots
#' @return List of plot objects
create_result3_figures <- function(result3_results, tree_name, output_dir) {
  cat("Creating Result 3 figures for", tree_name, "...\n")
  
  result <- result3_results$results[[tree_name]]
  
  plots <- list()
  
  # Plot trait signal for K
  if (!is.null(result$trait_analysis$signal_results)) {
    p_k <- plot_trait_signal(
      signal_results = result$trait_analysis$signal_results,
      metric_name = "K",
      title = paste("Blomberg's K by Subset Type -", tree_name)
    )
    
    plots$blombergs_k <- p_k
    
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, paste0("result3_blombergs_k_", tree_name, ".pdf"))
      ggsave(plot_file, p_k, width = 10, height = 6)
    }
    
    # Plot trait signal for Lambda
    p_lambda <- plot_trait_signal(
      signal_results = result$trait_analysis$signal_results,
      metric_name = "Lambda",
      title = paste("Pagel's λ by Subset Type -", tree_name)
    )
    
    plots$pagels_lambda <- p_lambda
    
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, paste0("result3_pagels_lambda_", tree_name, ".pdf"))
      ggsave(plot_file, p_lambda, width = 10, height = 6)
    }
  }
  
  # Plot pattern analysis
  if (!is.null(result$analysis_results$pattern_analysis)) {
    pattern_df <- result$analysis_results$pattern_analysis
    
    # Reshape for plotting
    plot_df <- data.frame(
      Model = rep(pattern_df$Model, 3),
      Subset_Type = rep(c("Clustered", "Random", "Dispersed"), each = nrow(pattern_df)),
      K_Mean = c(pattern_df$Clustered_Mean_K, pattern_df$Random_Mean_K, pattern_df$Dispersed_Mean_K),
      Lambda_Mean = c(pattern_df$Clustered_Mean_Lambda, pattern_df$Random_Mean_Lambda, pattern_df$Dispersed_Mean_Lambda)
    )
    
    # Plot K pattern
    p_k_pattern <- ggplot(plot_df, aes(x = Subset_Type, y = K_Mean, fill = Subset_Type)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      facet_wrap(~ Model, scales = "free_y") +
      geom_text(aes(label = round(K_Mean, 3)), vjust = -0.5, size = 3) +
      scale_fill_manual(values = c("Clustered" = "#CD5C5C", 
                                   "Random" = "#4682B4", 
                                   "Dispersed" = "#2E8B57")) +
      labs(title = paste("Mean Blomberg's K Pattern -", tree_name),
           x = "Subset Type", y = "Mean K") +
      theme_phylogenetic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots$k_pattern <- p_k_pattern
    
    # Plot Lambda pattern
    p_lambda_pattern <- ggplot(plot_df, aes(x = Subset_Type, y = Lambda_Mean, fill = Subset_Type)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      facet_wrap(~ Model, scales = "free_y") +
      geom_text(aes(label = round(Lambda_Mean, 3)), vjust = -0.5, size = 3) +
      scale_fill_manual(values = c("Clustered" = "#CD5C5C", 
                                   "Random" = "#4682B4", 
                                   "Dispersed" = "#2E8B57")) +
      labs(title = paste("Mean Pagel's λ Pattern -", tree_name),
           x = "Subset Type", y = "Mean λ") +
      theme_phylogenetic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots$lambda_pattern <- p_lambda_pattern
    
    # Save pattern plots
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, paste0("result3_k_pattern_", tree_name, ".pdf"))
      ggsave(plot_file, p_k_pattern, width = 10, height = 6)
      
      plot_file <- file.path(output_dir, paste0("result3_lambda_pattern_", tree_name, ".pdf"))
      ggsave(plot_file, p_lambda_pattern, width = 10, height = 6)
    }
  }
  
  return(plots)
}

#' Create multi-panel figure for Result 4
#'
#' @param result4_results Results from Result 4 analysis
#' @param output_dir Output directory for saving plots
#' @return List of plot objects
create_result4_figures <- function(result4_results, output_dir) {
  cat("Creating Result 4 figures...\n")
  
  result <- result4_results$result
  
  plots <- list()
  
  # Plot heuristic vs exact comparison
  if (!is.null(result$comparison$comparison_summary)) {
    p_comparison <- plot_heuristic_vs_exact(
      comparison_summary = result$comparison$comparison_summary,
      title = paste("Heuristic vs Exact Optimum -", result$tree_name)
    )
    
    plots$heuristic_vs_exact <- p_comparison
    
    # Plot gap analysis
    p_gap <- plot_gap_analysis(
      comparison_summary = result$comparison$comparison_summary,
      title = paste("Gap Analysis -", result$tree_name)
    )
    
    plots$gap_analysis <- p_gap
    
    # Save plots
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, "result4_heuristic_vs_exact.pdf")
      ggsave(plot_file, p_comparison, width = 8, height = 6)
      
      plot_file <- file.path(output_dir, "result4_gap_analysis.pdf")
      ggsave(plot_file, p_gap, width = 8, height = 6)
    }
  }
  
  # Plot rank analysis if available
  if (!is.null(result$evaluation$rank_analysis)) {
    rank_info <- result$evaluation$rank_analysis
    
    # Create rank visualization
    rank_df <- data.frame(
      Category = c("Better", "Equal", "Worse"),
      Count = c(rank_info$better_count, 
                rank_info$equal_count, 
                rank_info$evaluated_count - rank_info$better_count - rank_info$equal_count)
    )
    
    p_rank <- ggplot(rank_df, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_text(aes(label = Count), vjust = -0.5, size = 4) +
      scale_fill_manual(values = c("Better" = "#FF6B6B", 
                                   "Equal" = "#4ECDC4", 
                                   "Worse" = "#45B7D1")) +
      labs(title = paste("Rank Analysis -", result$tree_name),
           subtitle = paste("Percentile:", round(rank_info$percentile, 2), "%"),
           x = "Comparison", y = "Number of Subsets") +
      theme_phylogenetic()
    
    plots$rank_analysis <- p_rank
    
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, "result4_rank_analysis.pdf")
      ggsave(plot_file, p_rank, width = 8, height = 6)
    }
  }
  
  return(plots)
}

#' Create summary figure for all results
#'
#' @param all_results List containing results from all four analyses
#' @param cfg Configuration list
#' @param output_dir Output directory for saving plots
#' @return List of summary plot objects
create_summary_figures <- function(all_results, cfg, output_dir) {
  cat("Creating summary figures...\n")
  
  plots <- list()
  
  # Create output directory if it doesn't exist
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Result 1: Z-scores across trees
  if (!is.null(all_results$result1)) {
    result1 <- all_results$result1
    
    # Extract z-scores from all trees
    z_data <- data.frame()
    
    for (tree_name in names(result1$results)) {
      result <- result1$results[[tree_name]]
      z_scores <- result$summary$z_scores
      
      for (metric in names(z_scores)) {
        z_data <- rbind(z_data, data.frame(
          Tree = tree_name,
          Metric = metric,
          Z_Score = z_scores[[metric]],
          Significant = abs(z_scores[[metric]]) > 1.96,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(z_data) > 0) {
      p_z_summary <- ggplot(z_data, aes(x = Tree, y = Z_Score, fill = Significant)) +
        geom_bar(stat = "identity", alpha = 0.8) +
        facet_wrap(~ Metric, ncol = 3) +
        geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "red") +
        geom_text(aes(label = round(Z_Score, 2)), vjust = -0.5, size = 3) +
        scale_fill_manual(values = c("TRUE" = "#FF6B6B", "FALSE" = "#4ECDC4"),
                          name = "Significant (|Z| > 1.96)") +
        labs(title = "Result 1: Z-Scores Across Trees",
             x = "Tree", y = "Z-Score") +
        theme_phylogenetic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      plots$result1_z_summary <- p_z_summary
      
      if (!is.null(output_dir)) {
        plot_file <- file.path(output_dir, "summary_result1_z_scores.pdf")
        ggsave(plot_file, p_z_summary, width = 12, height = 6)
      }
    }
  }
  
  # Result 2: Algorithm performance across trees
  if (!is.null(all_results$result2)) {
    result2 <- all_results$result2
    
    # Extract algorithm metrics from all trees
    algo_data <- data.frame()
    
    for (tree_name in names(result2$results)) {
      result <- result2$results[[tree_name]]
      metrics <- result$comparison$metrics
      metrics$Tree <- tree_name
      algo_data <- rbind(algo_data, metrics)
    }
    
    if (nrow(algo_data) > 0) {
      # Plot for each metric
      for (metric in c("MinPD", "MeanPD", "MeanNND")) {
        p_algo <- ggplot(algo_data, aes(x = Tree, y = .data[[metric]], fill = Algorithm)) +
          geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
          labs(title = paste("Result 2:", metric, "by Algorithm Across Trees"),
               x = "Tree", y = metric) +
          theme_phylogenetic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        plots[[paste0("result2_", metric, "_summary")]] <- p_algo
        
        if (!is.null(output_dir)) {
          plot_file <- file.path(output_dir, paste0("summary_result2_", metric, ".pdf"))
          ggsave(plot_file, p_algo, width = 10, height = 6)
        }
      }
    }
  }
  
  # Result 3: Pattern consistency across trees
  if (!is.null(all_results$result3)) {
    result3 <- all_results$result3
    
    if (!is.null(result3$overall_summary) && nrow(result3$overall_summary) > 0) {
      # Plot pattern consistency
      pattern_summary <- result3$overall_summary %>%
        group_by(Model) %>%
        summarise(
          Pattern_K_Pct = mean(Pattern_K) * 100,
          Pattern_Lambda_Pct = mean(Pattern_Lambda) * 100,
          .groups = "drop"
        )
      
      # Reshape for plotting
      plot_df <- data.frame(
        Model = rep(pattern_summary$Model, 2),
        Metric = rep(c("K", "Lambda"), each = nrow(pattern_summary)),
        Pattern_Pct = c(pattern_summary$Pattern_K_Pct, pattern_summary$Pattern_Lambda_Pct)
      )
      
      p_pattern <- ggplot(plot_df, aes(x = Model, y = Pattern_Pct, fill = Metric)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
        geom_text(aes(label = paste0(round(Pattern_Pct, 1), "%")),
                  position = position_dodge(width = 0.8),
                  vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = c("K" = "#FF6B6B", "Lambda" = "#4ECDC4")) +
        labs(title = "Result 3: Pattern Consistency Across Trees",
             subtitle = "Percentage of trees where clustered > random > dispersed",
             x = "Model", y = "Pattern Holds (%)") +
        ylim(0, 100) +
        theme_phylogenetic()
      
      plots$result3_pattern_summary <- p_pattern
      
      if (!is.null(output_dir)) {
        plot_file <- file.path(output_dir, "summary_result3_pattern.pdf")
        ggsave(plot_file, p_pattern, width = 10, height = 6)
      }
    }
  }
  
  # Result 4: Heuristic performance
  if (!is.null(all_results$result4)) {
    result4 <- all_results$result4
    
    if (!is.null(result4$summary)) {
      # Create a performance summary
      perf_df <- data.frame(
        Metric = c("MinPD", "MeanPD", "MeanNND"),
        Rel_Gap = c(result4$summary$MinPD_Rel_Gap,
                    result4$summary$MeanPD_Rel_Gap,
                    result4$summary$MeanNND_Rel_Gap)
      )
      
      p_perf <- ggplot(perf_df, aes(x = Metric, y = abs(Rel_Gap), fill = Metric)) +
        geom_bar(stat = "identity", alpha = 0.8) +
        geom_text(aes(label = paste0(round(abs(Rel_Gap), 2), "%")),
                  vjust = -0.5, size = 4) +
        scale_fill_viridis(discrete = TRUE, option = "D") +
        labs(title = "Result 4: Heuristic Performance",
             subtitle = paste("Exact match:", result4$summary$Exact_Match),
             x = "Metric", y = "Absolute Relative Gap (%)") +
        theme_phylogenetic()
      
      plots$result4_performance_summary <- p_perf
      
      if (!is.null(output_dir)) {
        plot_file <- file.path(output_dir, "summary_result4_performance.pdf")
        ggsave(plot_file, p_perf, width = 8, height = 6)
      }
    }
  }
  
  return(plots)
}

#' Test plotting functions
test_plotting <- function() {
  # Load required libraries
  library(ggplot2)
  library(gridExtra)
  library(viridis)
  library(dplyr)
  
  cat("Testing plotting functions...\n")
  
  # Create test data
  set.seed(123)
  
  # Test null distribution plot
  cat("1. Testing null distribution plot...\n")
  null_data <- data.frame(
    MinPD = rnorm(1000, mean = 10, sd = 2),
    MeanPD = rnorm(1000, mean = 20, sd = 3),
    MeanNND = rnorm(1000, mean = 5, sd = 1)
  )
  
  observed <- list(MinPD = 15, MeanPD = 25, MeanNND = 7)
  
  p_null <- plot_null_distribution(null_data, observed, "MinPD", 
                                   title = "Test Null Distribution")
  print(p_null)
  
  # Test algorithm comparison plot
  cat("2. Testing algorithm comparison plot...\n")
  algo_data <- data.frame(
    Algorithm = c("A_main", "B_greedy", "C_meanpd", "D_minpd", "E_meannnd"),
    MinPD = c(15, 14, 12, 16, 13),
    MeanPD = c(25, 23, 26, 22, 24),
    MeanNND = c(7, 6, 5, 8, 6)
  )
  
  algo_plots <- plot_algorithm_comparison(algo_data, title = "Test Algorithm Comparison")
  print(algo_plots$combined)
  
  # Test trait signal plot
  cat("3. Testing trait signal plot...\n")
  signal_data <- data.frame(
    Model = rep(c("BM", "OU_alpha0.2", "OU_alpha1"), each = 30),
    Subset_Type = rep(rep(c("dispersed", "clustered", "random"), each = 10), 3),
    K = c(rnorm(10, 0.8, 0.1), rnorm(10, 1.2, 0.1), rnorm(10, 1.0, 0.1),
          rnorm(10, 0.7, 0.1), rnorm(10, 1.1, 0.1), rnorm(10, 0.9, 0.1),
          rnorm(10, 0.6, 0.1), rnorm(10, 1.0, 0.1), rnorm(10, 0.8, 0.1)),
    Lambda = c(rnorm(10, 0.7, 0.1), rnorm(10, 1.1, 0.1), rnorm(10, 0.9, 0.1),
               rnorm(10, 0.6, 0.1), rnorm(10, 1.0, 0.1), rnorm(10, 0.8, 0.1),
               rnorm(10, 0.5, 0.1), rnorm(10, 0.9, 0.1), rnorm(10, 0.7, 0.1))
  )
  
  p_signal <- plot_trait_signal(signal_data, metric_name = "K", 
                                title = "Test Trait Signal")
  print(p_signal)
  
  # Test heuristic vs exact plot
  cat("4. Testing heuristic vs exact plot...\n")
  comp_data <- data.frame(
    Metric = c("MinPD", "MeanPD", "MeanNND"),
    Heuristic = c(15, 25, 7),
    Exact = c(16, 26, 8),
    Gap = c(1, 1, 1),
    Rel_Gap_Pct = c(6.25, 3.85, 12.5)
  )
  
  p_heuristic <- plot_heuristic_vs_exact(comp_data, title = "Test Heuristic vs Exact")
  print(p_heuristic)
  
  p_gap <- plot_gap_analysis(comp_data, title = "Test Gap Analysis")
  print(p_gap)
  
  cat("\nAll plotting tests completed successfully.\n")
  
  return(list(
    null_plot = p_null,
    algo_plots = algo_plots,
    signal_plot = p_signal,
    heuristic_plot = p_heuristic,
    gap_plot = p_gap
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_plotting()
  cat("\nPlotting test completed successfully.\n")
}
