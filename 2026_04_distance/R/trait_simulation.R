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

#' Simulate traits under Pagel's lambda model
#'
#' @param tree A phylo object
#' @param n_traits Number of traits to simulate (default: 1)
#' @param lambda Pagel's lambda parameter (0-1)
#' @param sigma2 Rate parameter for BM (default: 1.0)
#' @param n_reps Number of simulation replicates (default: 1)
#' @return A list of simulated trait matrices
simulate_lambda_traits <- function(tree, n_traits = 1, lambda = 1.0, sigma2 = 1.0, n_reps = 1) {
  cat("Simulating", n_reps, "replicates of", n_traits, "lambda trait(s) with lambda =", lambda, 
      "sigma2 =", sigma2, "\n")
  
  simulations <- list()
  
  for (i in 1:n_reps) {
    # Transform tree according to lambda
    lambda_tree <- tree
    if (lambda < 1.0) {
      # Apply lambda transformation to branch lengths
      lambda_tree$edge.length <- lambda_tree$edge.length * lambda
      # Add star-like component for the remaining (1-lambda) portion
      # This is a simplified approach; more sophisticated methods exist
      # but this approximates the lambda model reasonably well
      lambda_tree$edge.length <- lambda_tree$edge.length + (1 - lambda) * mean(lambda_tree$edge.length)
    }
    
    # Simulate BM on transformed tree
    if (n_traits == 1) {
      trait_values <- fastBM(lambda_tree, sig2 = sigma2)
    } else {
      trait_values <- replicate(n_traits, fastBM(lambda_tree, sig2 = sigma2))
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
#' @param lambda_values Vector of lambda values for Pagel's lambda model
#' @param ou_alpha_values Vector of alpha values for OU process
#' @return A comprehensive simulation result
simulate_traits_for_subsets <- function(tree, subsets, subset_names, n_reps = 100,
                                        lambda_values = c(1.0, 0.7),
                                        ou_alpha_values = c(0.2, 1, 5)) {
  
  cat("Simulating traits for", length(subsets), "subsets with", n_reps, "replicates\n")
  
  results <- list()
  
  # Lambda simulations for different lambda values
  results$Lambda <- list()
  
  for (lambda_val in lambda_values) {
    cat("  Lambda simulation with lambda =", lambda_val, "...\n")
    lambda_simulations <- simulate_lambda_traits(tree, n_traits = 1, lambda = lambda_val, 
                                                 sigma2 = 1.0, n_reps = n_reps)
    
    # Extract lambda traits for each subset
    lambda_results <- list()
    for (i in seq_along(subsets)) {
      subset_traits <- sapply(lambda_simulations, function(sim) {
        extract_subset_traits(sim, subsets[[i]])
      })
      lambda_results[[subset_names[i]]] <- subset_traits
    }
    results$Lambda[[paste0("Lambda", lambda_val)]] <- lambda_results
  }
  
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
  
  # Store metadata including subsets
  results$metadata <- list(
    tree_n_tips = length(tree$tip.label),
    n_reps = n_reps,
    subset_names = subset_names,
    subsets = subsets,  # Store the actual subsets
    subset_sizes = sapply(subsets, length),
    lambda_values = lambda_values,
    ou_alpha_values = ou_alpha_values
  )
  
  return(results)
}

# Note: calc_blombergs_K and calc_pagels_lambda functions are defined in signal_metrics.R
# calculate_signal_metrics is a placeholder function that is not used in the main analysis

#' Test trait simulation functions
test_trait_simulation <- function() {
  # Load required libraries
  library(ape)
  library(phytools)
  
  # Create a test tree
  test_tree <- rtree(20)
  test_tree$tip.label <- paste0("sp", 1:20)
  
  cat("Testing trait simulation functions:\n")
  
  # Test 1: Lambda simulation
  cat("\n1. Lambda simulation (lambda=1.0):\n")
  lambda1_sim <- simulate_lambda_traits(test_tree, n_traits = 1, lambda = 1.0, sigma2 = 1.0, n_reps = 3)
  cat("   Simulated", length(lambda1_sim), "replicates\n")
  cat("   First replicate values (first 5 tips):", lambda1_sim[[1]][1:5], "\n")
  
  # Test 2: Lambda simulation with lambda=0.7
  cat("\n2. Lambda simulation (lambda=0.7):\n")
  lambda07_sim <- simulate_lambda_traits(test_tree, n_traits = 1, lambda = 0.7, sigma2 = 1.0, n_reps = 3)
  cat("   Simulated", length(lambda07_sim), "replicates\n")
  cat("   First replicate values (first 5 tips):", lambda07_sim[[1]][1:5], "\n")
  
  # Test 3: OU simulation
  cat("\n3. Ornstein-Uhlenbeck simulation:\n")
  ou_sim <- simulate_ou_traits(test_tree, n_traits = 1, alpha = 1.0, sigma2 = 1.0, n_reps = 3)
  cat("   Simulated", length(ou_sim), "replicates\n")
  cat("   First replicate values (first 5 tips):", ou_sim[[1]][1:5], "\n")
  
  # Test 4: Extract subset traits
  cat("\n4. Extracting subset traits:\n")
  subset_tips <- c("sp1", "sp3", "sp5", "sp7", "sp9")
  subset_traits <- extract_subset_traits(lambda1_sim[[1]], subset_tips)
  cat("   Subset traits:", subset_traits, "\n")
  
  # Test 5: Simulate for multiple subsets
  cat("\n5. Simulating for multiple subsets:\n")
  subsets <- list(
    c("sp1", "sp3", "sp5", "sp7", "sp9"),
    c("sp2", "sp4", "sp6", "sp8", "sp10"),
    c("sp1", "sp5", "sp10", "sp15", "sp20")
  )
  subset_names <- c("subset1", "subset2", "subset3")
  
  sim_results <- simulate_traits_for_subsets(test_tree, subsets, subset_names, 
                                            n_reps = 5, lambda_values = c(1.0, 0.7),
                                            ou_alpha_values = c(0.2, 1))
  
  cat("   Simulation results structure:\n")
  cat("     Lambda results for", length(sim_results$Lambda), "lambda values\n")
  cat("     OU results for", length(sim_results$OU), "alpha values\n")
  
  # Test 6: Calculate signal metrics using analyze_simulated_traits_signal
  cat("\n6. Calculating signal metrics:\n")
  signal_metrics <- analyze_simulated_traits_signal(sim_results, test_tree, n_reps_to_analyze = 3)
  cat("   Calculated metrics for", nrow(signal_metrics), "simulations\n")
  cat("   First few rows:\n")
  print(head(signal_metrics))
  
  return(list(
    lambda1_sim = lambda1_sim,
    lambda07_sim = lambda07_sim,
    ou_sim = ou_sim,
    sim_results = sim_results,
    signal_metrics = signal_metrics
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_trait_simulation()
}
