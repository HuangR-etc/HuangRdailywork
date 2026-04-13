# Phylogenetic signal metrics calculation
# This module calculates Blomberg's K and Pagel's lambda for trait data

library(ape)
library(picante)
library(phytools)

#' Calculate Blomberg's K for trait data
#'
#' @param trait_values Vector or matrix of trait values
#' @param tree A phylo object
#' @return Blomberg's K value
calc_blombergs_K <- function(trait_values, tree) {
  # Ensure tree tips and trait names match
  if (is.vector(trait_values)) {
    # Single trait
    trait_values <- trait_values[tree$tip.label]
    
    # Calculate K using picante
    K_result <- Kcalc(trait_values, tree)
    
    return(K_result)
  } else {
    # Multiple traits
    n_traits <- ncol(trait_values)
    K_values <- numeric(n_traits)
    
    for (i in 1:n_traits) {
      trait_vec <- trait_values[tree$tip.label, i]
      K_values[i] <- Kcalc(trait_vec, tree)
    }
    
    names(K_values) <- colnames(trait_values)
    return(K_values)
  }
}

#' Calculate Pagel's lambda for trait data
#'
#' @param trait_values Vector or matrix of trait values
#' @param tree A phylo object
#' @return Pagel's lambda value
calc_pagels_lambda <- function(trait_values, tree) {
  # Ensure tree tips and trait names match
  if (is.vector(trait_values)) {
    # Single trait
    trait_values <- trait_values[tree$tip.label]
    
    # Calculate lambda using phytools
    lambda_result <- phylosig(tree, trait_values, method = "lambda", test = FALSE)
    
    # phylosig returns a list with lambda value when method = "lambda"
    return(lambda_result$lambda)
  } else {
    # Multiple traits
    n_traits <- ncol(trait_values)
    lambda_values <- numeric(n_traits)
    
    for (i in 1:n_traits) {
      trait_vec <- trait_values[tree$tip.label, i]
      lambda_result <- phylosig(tree, trait_vec, method = "lambda", test = FALSE)
      lambda_values[i] <- lambda_result$lambda
    }
    
    names(lambda_values) <- colnames(trait_values)
    return(lambda_values)
  }
}

#' Calculate Moran's I for trait data with phylogenetic weights
#'
#' @param trait_values Vector of trait values
#' @param tree A phylo object
#' @return A list with Moran's I statistics: observed, expected, sd, p.value
calc_morans_I <- function(trait_values, tree) {
  # Ensure tree tips and trait names match
  trait_values <- trait_values[tree$tip.label]
  
  # Calculate patristic distance matrix
  dist_matrix <- cophenetic.phylo(tree)
  
  # Construct weight matrix W: W[i,j] = 1/D[i,j] when i != j, W[i,i] = 0
  W <- 1 / dist_matrix
  diag(W) <- 0
  
  # Ensure W is symmetric and finite
  W[!is.finite(W)] <- 0
  
  # Calculate Moran's I using ape::Moran.I
  moran_result <- ape::Moran.I(trait_values, weight = W)
  
  return(moran_result)
}

#' Calculate phylogenetic signal for multiple subsets
#'
#' @param tree Full tree
#' @param subsets List of subsets (each is vector of tip names)
#' @param subset_names Names for the subsets
#' @param trait_values Trait values for all tips
#' @param metrics Vector of metrics to calculate ("K", "lambda", or both)
#' @return A data frame with signal metrics for all subsets
calculate_subsets_signal <- function(tree, subsets, subset_names, trait_values, 
                                     metrics = c("K", "lambda")) {
  
  cat("Calculating phylogenetic signal for", length(subsets), "subsets\n")
  
  results <- data.frame()
  
  for (i in seq_along(subsets)) {
    subset_tips <- subsets[[i]]
    subset_name <- subset_names[i]
    
    cat("  Processing subset:", subset_name, "(", length(subset_tips), "tips)\n")
    
    # Extract subtree for this subset
    subset_tree <- keep.tip(tree, subset_tips)
    
    # Extract trait values for this subset
    if (is.vector(trait_values)) {
      subset_trait_values <- trait_values[subset_tips]
    } else {
      subset_trait_values <- trait_values[subset_tips, , drop = FALSE]
    }
    
    # Calculate requested metrics
    subset_result <- data.frame(
      Subset = subset_name,
      N_Tips = length(subset_tips),
      stringsAsFactors = FALSE
    )
    
    if ("K" %in% metrics) {
      tryCatch({
        K_value <- calc_blombergs_K(subset_trait_values, subset_tree)
        if (is.vector(K_value) && length(K_value) == 1) {
          subset_result$K <- K_value
        } else {
          # For multiple traits, take the mean
          subset_result$K <- mean(K_value, na.rm = TRUE)
        }
      }, error = function(e) {
        cat("    Error calculating K:", e$message, "\n")
        subset_result$K <- NA
      })
    }
    
    if ("lambda" %in% metrics) {
      tryCatch({
        lambda_value <- calc_pagels_lambda(subset_trait_values, subset_tree)
        if (is.vector(lambda_value) && length(lambda_value) == 1) {
          subset_result$Lambda <- lambda_value
        } else {
          # For multiple traits, take the mean
          subset_result$Lambda <- mean(lambda_value, na.rm = TRUE)
        }
      }, error = function(e) {
        cat("    Error calculating lambda:", e$message, "\n")
        subset_result$Lambda <- NA
      })
    }
    
    results <- rbind(results, subset_result)
  }
  
  return(results)
}

#' Calculate signal metrics for simulated traits
#'
#' @param simulation_results Results from simulate_traits_for_subsets
#' @param tree Full tree
#' @param n_reps_to_analyze Number of replicates to analyze (for performance)
#' @return A comprehensive data frame with signal metrics
analyze_simulated_traits_signal <- function(simulation_results, tree, n_reps_to_analyze = 100) {
  cat("Analyzing phylogenetic signal in simulated traits...\n")
  
  all_results <- data.frame()
  
  # Get subset information from metadata
  subset_names <- simulation_results$metadata$subset_names
  subsets <- simulation_results$metadata$subsets
  
  if (is.null(subsets)) {
    warning("No subsets found in simulation results metadata. Cannot calculate signal metrics.")
    return(all_results)
  }
  
  # Process Lambda simulations
  lambda_models <- names(simulation_results$Lambda)
  
  for (lambda_model in lambda_models) {
    lambda_results <- simulation_results$Lambda[[lambda_model]]
    
    for (subset_name in names(lambda_results)) {
      cat("  Processing", lambda_model, "for subset:", subset_name, "\n")
      
      # Find the index of this subset
      subset_idx <- which(subset_names == subset_name)
      if (length(subset_idx) == 0) {
        warning(paste("Subset", subset_name, "not found in metadata"))
        next
      }
      
      # Get the tip names for this subset
      subset_tips <- subsets[[subset_idx]]
      
      # Get trait values for this subset
      subset_traits <- lambda_results[[subset_name]]
      n_reps <- min(n_reps_to_analyze, ncol(subset_traits))
      
      for (rep in 1:n_reps) {
        # Extract trait values for this replicate
        trait_values <- subset_traits[, rep]
        names(trait_values) <- subset_tips
        
        # Extract subtree for this subset
        subset_tree <- keep.tip(tree, subset_tips)
        
        # Calculate K
        K_value <- tryCatch({
          calc_blombergs_K(trait_values, subset_tree)
        }, error = function(e) {
          cat("    Error calculating K for replicate", rep, ":", e$message, "\n")
          NA
        })
        
        # Calculate Lambda
        lambda_value <- tryCatch({
          calc_pagels_lambda(trait_values, subset_tree)
        }, error = function(e) {
          cat("    Error calculating Lambda for replicate", rep, ":", e$message, "\n")
          NA
        })
        
        # Calculate Moran's I
        moran_result <- tryCatch({
          calc_morans_I(trait_values, subset_tree)
        }, error = function(e) {
          cat("    Error calculating Moran's I for replicate", rep, ":", e$message, "\n")
          list(observed = NA, expected = NA, sd = NA, p.value = NA)
        })
        
        all_results <- rbind(all_results, data.frame(
          Model = lambda_model,
          Subset = subset_name,
          Replicate = rep,
          K = K_value,
          Lambda = lambda_value,
          MoransI = moran_result$observed,
          MoransI_expected = moran_result$expected,
          MoransI_sd = moran_result$sd,
          MoransI_p = moran_result$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Process OU simulations
  ou_models <- names(simulation_results$OU)
  
  for (ou_model in ou_models) {
    ou_results <- simulation_results$OU[[ou_model]]
    
    for (subset_name in names(ou_results)) {
      cat("  Processing", ou_model, "for subset:", subset_name, "\n")
      
      # Find the index of this subset
      subset_idx <- which(subset_names == subset_name)
      if (length(subset_idx) == 0) {
        warning(paste("Subset", subset_name, "not found in metadata"))
        next
      }
      
      # Get the tip names for this subset
      subset_tips <- subsets[[subset_idx]]
      
      # Get trait values for this subset
      subset_traits <- ou_results[[subset_name]]
      n_reps <- min(n_reps_to_analyze, ncol(subset_traits))
      
      for (rep in 1:n_reps) {
        # Extract trait values for this replicate
        trait_values <- subset_traits[, rep]
        names(trait_values) <- subset_tips
        
        # Extract subtree for this subset
        subset_tree <- keep.tip(tree, subset_tips)
        
        # Calculate K
        K_value <- tryCatch({
          calc_blombergs_K(trait_values, subset_tree)
        }, error = function(e) {
          cat("    Error calculating K for replicate", rep, ":", e$message, "\n")
          NA
        })
        
        # Calculate Lambda
        lambda_value <- tryCatch({
          calc_pagels_lambda(trait_values, subset_tree)
        }, error = function(e) {
          cat("    Error calculating Lambda for replicate", rep, ":", e$message, "\n")
          NA
        })
        
        # Calculate Moran's I
        moran_result <- tryCatch({
          calc_morans_I(trait_values, subset_tree)
        }, error = function(e) {
          cat("    Error calculating Moran's I for replicate", rep, ":", e$message, "\n")
          list(observed = NA, expected = NA, sd = NA, p.value = NA)
        })
        
        all_results <- rbind(all_results, data.frame(
          Model = ou_model,
          Subset = subset_name,
          Replicate = rep,
          K = K_value,
          Lambda = lambda_value,
          MoransI = moran_result$observed,
          MoransI_expected = moran_result$expected,
          MoransI_sd = moran_result$sd,
          MoransI_p = moran_result$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(all_results)
}

#' Test signal metrics functions
test_signal_metrics <- function() {
  # Load required libraries
  library(ape)
  library(picante)
  library(phytools)
  
  # Create a test tree
  test_tree <- rtree(20)
  test_tree$tip.label <- paste0("sp", 1:20)
  
  # Simulate trait data under BM
  set.seed(123)
  trait_values <- fastBM(test_tree, sig2 = 1.0)
  
  cat("Testing signal metrics functions:\n")
  
  # Test 1: Calculate Blomberg's K for full tree
  cat("\n1. Calculating Blomberg's K for full tree:\n")
  K_full <- calc_blombergs_K(trait_values, test_tree)
  cat("   K =", K_full, "\n")
  
  # Test 2: Calculate Pagel's lambda for full tree
  cat("\n2. Calculating Pagel's lambda for full tree:\n")
  lambda_full <- calc_pagels_lambda(trait_values, test_tree)
  cat("   Lambda =", lambda_full, "\n")
  
  # Test 3: Calculate signal for subsets
  cat("\n3. Calculating signal for subsets:\n")
  subsets <- list(
    c("sp1", "sp3", "sp5", "sp7", "sp9"),
    c("sp2", "sp4", "sp6", "sp8", "sp10"),
    c("sp1", "sp5", "sp10", "sp15", "sp20")
  )
  subset_names <- c("subset1", "subset2", "subset3")
  
  subset_signal <- calculate_subsets_signal(test_tree, subsets, subset_names, 
                                           trait_values, metrics = c("K", "lambda"))
  
  cat("   Subset signal metrics:\n")
  print(subset_signal)
  
  # Test 4: Simulate traits and analyze signal
  cat("\n4. Simulating traits and analyzing signal (simplified):\n")
  
  # Create a simple simulation result structure
  sim_results <- list(
    BM = list(
      subset1 = matrix(rnorm(5*10), nrow = 5, ncol = 10),
      subset2 = matrix(rnorm(5*10), nrow = 5, ncol = 10),
      subset3 = matrix(rnorm(5*10), nrow = 5, ncol = 10)
    ),
    OU = list(
      OU_alpha0.2 = list(
        subset1 = matrix(rnorm(5*10), nrow = 5, ncol = 10),
        subset2 = matrix(rnorm(5*10), nrow = 5, ncol = 10),
        subset3 = matrix(rnorm(5*10), nrow = 5, ncol = 10)
      )
    ),
    metadata = list(
      subset_names = c("subset1", "subset2", "subset3"),
      subsets = subsets
    )
  )
  
  signal_analysis <- analyze_simulated_traits_signal(sim_results, test_tree, n_reps_to_analyze = 5)
  
  cat("   Signal analysis results (first few rows):\n")
  print(head(signal_analysis))
  
  return(list(
    K_full = K_full,
    lambda_full = lambda_full,
    subset_signal = subset_signal,
    signal_analysis = signal_analysis
  ))
}

# Run test if executed directly
if (sys.nframe() == 0) {
  test_results <- test_signal_metrics()
}
