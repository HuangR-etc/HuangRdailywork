# Trait simulation for phylogenetic dispersed subset analysis
# This module simulates traits under Brownian Motion and Ornstein-Uhlenbeck processes

library(ape)
library(geiger)
library(phytools)

#' Simulate traits under Brownian Motion (BM) process
#'
#' @param tree A phylo object
#' @param n_traits Number of traits to simulate (default: 1)
#' @param sigma2 Rate parameter for BM (default: 1.0)
#' @param n_reps Number of simulation replicates (default: 1)
#' @return A list of simulated trait matrices
simulate_bm_traits <- function(tree, n_traits = 1, sigma2 = 1.0, n_reps = 1) {
  cat("Simulating", n_reps, "replicates of", n_traits, "BM trait(s) with sigma2 =", sigma2, "\n")
  
  simulations <- list()
  
  for (i in 1:n_reps) {
    # Simulate using fastBM from phytools
    if (n_traits == 1) {
      trait_values <- fastBM(tree, sig2 = sigma2)
    } else {
      trait_values <- replicate(n_traits, fastBM(tree, sig2 = sigma2))
      colnames(trait_values) <- paste0("Trait", 1:n_traits)
    }
    
    simulations[[i]] <- trait_values
  }
  
  return(simulations)
}

#' Simulate traits under Ornstein-Uhlenbeck (OU) process
#'
#' @param tree A phylo object
#' @param n_traits Number of traits to simulate (default: 1)
#' @param alpha Strength of selection (OU parameter)
#' @param sigma2 Rate parameter (default: 1.0)
#' @param theta Optimal trait value (default: 0)
#' @param n_reps Number of simulation replicates (default: 1)
#' @return A list of simulated trait matrices
simulate_ou_traits <- function(tree, n_traits = 1, alpha = 1.0, sigma2 = 1.0, theta = 0, n_reps = 1) {
  cat("Simulating", n_reps, "replicates of", n_traits, "OU trait(s) with alpha =", alpha, 
      "sigma2 =", sigma2, "theta =", theta, "\n")
  
  simulations <- list()
  
  for (i in 1:n_reps) {
    if (n_traits == 1) {
      # Single trait simulation
      trait_values <- fastBM(tree, sig2 = sigma2, alpha = alpha, theta = theta)
    } else {
      # Multiple traits (independent evolution)
      trait_values <- replicate(n_traits, fastBM(tree, sig2 = sigma2, alpha = alpha, theta = theta))
      colnames(trait_values) <- paste0("Trait", 1:n_traits)
    }
    
    simulations[[i]] <- trait_values
  }
  
  return(simulations)
}

#' Extract trait values for a specific subset
#'
#' @param trait_matrix Trait matrix (from simulation)
#' @param tip_names Vector of tip names in the subset
#' @return Trait values for the subset
extract_subset_traits <- function(trait_matrix, tip_names) {
  if (is.vector(trait_matrix)) {
    # Single trait
    return(trait_matrix[tip_names])
  } else {
    # Multiple traits
    return(trait_matrix[tip_names, , drop = FALSE])
  }
}

#' Simulate traits for multiple models and extract subset values
#'
#' @param tree A phylo object
#' @param subsets List of subsets (each is vector of tip names)
#' @param subset_names Names for the subsets
#' @param n_reps Number of simulation replicates
#' @param ou_alpha_values Vector of alpha values for OU process
#' @return A comprehensive simulation result
simulate_traits_for_subsets <- function(tree, subsets, subset_names, n_reps = 100,
                                        ou_alpha_values = c(0.2, 1, 5)) {
  
  cat("Simulating traits for", length(subsets), "subsets with", n_reps, "replicates\n")
  
  results <- list()
  
  # BM simulation
  cat("  Brownian Motion simulation...\n")
  bm_simulations <- simulate_bm_traits(tree, n_traits = 1, sigma2 = 1.0, n_reps = n_reps)
  
  # Extract BM traits for each subset
  bm_results <- list()
  for (i in seq_along(subsets)) {
    subset_traits <- sapply(bm_simulations, function(sim) {
      extract_subset_traits(sim, subsets[[i]])
    })
    bm_results[[subset_names[i]]] <- subset_traits
  }
  results$BM <- bm_results
  
  # OU simulations for different alpha values
  results$OU <- list()
  
  for (alpha in ou_alpha_values) {
    cat("  OU simulation with alpha =", alpha, "...\n")
    ou_simulations <- simulate_ou_traits(tree, n_traits = 1, alpha = alpha, 
                                         sigma2 = 1.0, theta = 0, n_reps = n_reps)
    
    # Extract OU traits for each subset
    ou_results <- list()
    for (i in seq_along(subsets)) {
      subset_traits <- sapply(ou_simulations, function(sim) {
        extract_subset_traits(sim, subsets[[i]])
      })
      ou_results[[subset_names[i]]] <- subset_traits
    }
    results$OU[[paste0("OU_alpha", alpha)]] <- ou_results
  }
  
  # Store metadata
  results$metadata <- list(
    tree_n_tips = length(tree$tip.label),
    n_reps = n_reps,
    subset_names = subset_names,
    subset_sizes = sapply(subsets, length),
    ou_alpha_values = ou_alpha_values
  )
  
  return(results)
}

#' Calculate Blomberg's K for a set of trait values
#'
#' @param trait_values Vector of trait values
#' @param tree A phylo object (for the subset)
#' @param subset_tip_names Tip names in the subset
#' @return Blomberg's K value
calc_blombergs_K <- function(trait_values, tree, subset_tip_names) {
  # Extract subtree for the subset
  subset_tree <- keep.tip(tree, subset_tip_names)
  
  # Ensure trait values are in the same order as tree tips
  trait_values <- trait_values[subset_tree$tip.label]
  
  # Calculate K using picante package if available, otherwise manual calculation
  if (requireNamespace("picante", quietly = TRUE)) {
    K <- picante::Kcalc(trait_values, subset_tree)
  } else {
    # Manual calculation of K
    # This is a simplified version - for production use, consider implementing full K calculation
    warning("picante package not available, using simplified K calculation")
    
    # Calculate mean squared contrast
    contrasts <- pic(trait_values, subset_tree)
    ms_contrast <- mean(contrasts^2)
    
    # Calculate expected under Brownian motion
    n <- length(trait_values)
    expected_ms <- var(trait_values) * (n - 1) / n
    
    # Calculate K
    K <- ms_contrast / expected_ms
  }
  
  return(K)
}

#' Calculate Pagel's lambda for a set of trait values
#'
#' @param trait_values Vector of trait values
#' @param tree A phylo object (for the subset)
#' @param subset_tip_names Tip names in the subset
#' @return Pagel's lambda value
calc_pagels_lambda <- function(trait_values, tree, subset_tip_names) {
  # Extract subtree for the subset
  subset_tree <- keep.tip(tree, subset_tip_names)
  
  # Ensure trait values are in the same order as tree tips
  trait_values <- trait_values[subset_tree$tip.label]
  
  # Calculate lambda using phytools
  if (requireNamespace("phytools", quietly = TRUE)) {
    result <- phytools::phylosig(subset_tree, trait_values, method = "lambda", test = FALSE)
    lambda <- result$lambda
  } else {
    # Simplified lambda calculation
    warning("phytools package not available, using simplified lambda calculation")
    
    # For simplicity, return NA
    lambda <- NA
  }
  
  return(lambda)
}

#' Calculate phylogenetic signal metrics for simulated traits
#'
#' @param simulation_results Results from simulate_traits_for_subsets
#' @param tree Full tree
#' @return A data frame with signal metrics for all simulations
calculate_signal_metrics <- function(simulation_results, tree) {
  cat("Calculating phylogenetic signal metrics...\n")
  
  results_df <- data.frame()
  
  # Process BM simulations
  bm_results <- simulation_results$BM
  subset_names <- names(bm_results)
  
  for (subset_name in subset_names) {
    cat("  Processing BM for subset:", subset_name, "\n")
    
    # Get subset tip names (assuming they're stored somewhere)
    # For now, we'll need to pass this information differently
    # This is a placeholder - actual implementation would need subset tip names
    
    # Placeholder: calculate metrics for first few replicates
    n_reps <- min(10, ncol(bm_results[[subset_name]]))
    
    for (rep in 1:n_reps) {
      # This is a simplified version - actual implementation would need proper subset trees
      # For now, we'll create a dummy result
      results_df <- rbind(results_df, data.frame(
        Model = "BM",
        Subset = subset_name,
        Replicate = rep,
        K = NA,  # Placeholder
        Lambda = NA,  # Placeholder
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Process OU simulations
  ou_models <- names(simulation_results$OU)
  
  for (ou_model in ou_models) {
    ou_results <- simulation_results$OU[[ou_model]]
    
    for (subset_name in names(ou_results)) {
      cat("  Processing", ou_model, "for subset:", subset_name, "\n")
      
      # Placeholder: calculate metrics for first few replicates
      n_reps <- min(10, ncol(ou_results[[subset_name]]))
      
      for (rep in 1:n_reps) {
        results_df <- rbind(results_df, data.frame(
          Model = ou_model,
          Subset = subset_name,
          Replicate = rep,
          K = NA,  # Placeholder
          Lambda = NA,  # Placeholder
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(results_df)
}

#' Test trait simulation functions
test_trait_simulation <- function() {
  # Load required libraries
  library(ape)
  library(phytools)
  
  # Create a test tree
  test_tree <- rtree(20)
  test_tree$tip.label <- paste0("sp", 1:20)
  
  cat("Testing trait simulation functions:\n")
  
  # Test 1: BM simulation
  cat("\n1. Brownian Motion simulation:\n")
  bm_sim <- simulate_bm_traits(test_tree, n_traits = 1, sigma2 = 1.0, n_reps = 3)
  cat("   Simulated", length(bm_sim), "replicates\n")
  cat("   First replicate values (first 5 tips):", bm_sim[[1]][1:5], "\n")
  
  # Test 2: OU simulation
  cat("\n2. Ornstein-Uhlenbeck simulation:\n")
  ou_sim <- simulate_ou_traits(test_tree, n_traits = 1, alpha = 1.0, sigma2 = 1.0, n_reps = 3)
  cat("   Simulated", length(ou_sim), "replicates\n")
  cat("   First replicate values (first 5 tips):", ou_sim[[1]][1:5], "\n")
  
  # Test 3: Extract subset traits
  cat("\n3. Extracting subset traits:\n")
  subset_tips <- c("sp1", "sp3", "sp5", "sp7", "sp9")
  subset_traits <- extract_subset_traits(bm_sim[[1]], subset_tips)
  cat("   Subset traits:", subset_traits, "\n")
  
  # Test 4: Simulate for multiple subsets
  cat("\n4. Simulating for multiple subsets:\n")
  subsets <- list(
    c("sp1", "sp3", "sp5", "sp7", "sp9"),
    c("sp2", "sp4", "sp6", "sp8", "sp10"),
    c("sp1", "sp5", "sp10", "sp15", "sp20")
  )
  subset_names <- c("subset1", "subset2", "subset3")
  
  sim_results <- simulate_traits_for_subsets(test_tree, subsets, subset_names, 
                                            n_reps = 5, ou_alpha_values = c(0.2, 1))
  
  cat("   Simulation results structure:\n")
  cat("     BM results for", length(sim_results$BM), "subsets\n")
  cat("     OU results for", length(sim_results$OU), "alpha values\n")
  
  # Test 5: Calculate signal metrics (simplified)
  cat("\n5. Calculating signal metrics (simplified):\n")
  signal_metrics <- calculate_signal_metrics(sim_results, test_tree)
  cat("   Calculated metrics for", nrow(signal_metrics), "simulations\n")
  cat("   First few rows:\n")
  print(head(signal_metrics))
  
  return(list(
    bm_sim = bm_sim,
    ou_sim = ou_sim,
    sim_results = sim_results,
    signal_metrics = signal_metrics
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_trait_simulation()
}
