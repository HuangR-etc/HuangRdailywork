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

#' Calculate Mean Phylogenetic Nearest-Neighbor Similarity (mpnns)
#'
#' @param species Vector of species names in the subset
#' @param trait_values Vector of trait values aligned with species
#' @param phylo_dist_matrix Phylogenetic distance matrix for all species
#' @return mpnns value (NA if subset size < 2 or insufficient valid data)
calc_mpnns <- function(species, trait_values, phylo_dist_matrix) {
  # species: current subset species name vector
  # trait_values: numeric vector aligned with species
  # phylo_dist_matrix: phylogenetic distance matrix for all species
  
  # Check if subset size is less than 2
  if (length(species) < 2) {
    return(NA)
  }
  
  # 1. Extract distance submatrix for the subset
  dsub <- phylo_dist_matrix[species, species]
  
  # 2. Set diagonal to Inf to avoid self as nearest neighbor
  diag(dsub) <- Inf
  
  # 3. For each species, find the nearest neighbor index
  nn_idx <- apply(dsub, 1, which.min)
  
  # 4. Calculate similarity for each species with its nearest neighbor
  sims <- numeric(length(species))
  for (i in seq_along(species)) {
    j <- nn_idx[i]
    # Check for NA trait values
    if (is.na(trait_values[i]) || is.na(trait_values[j])) {
      sims[i] <- NA
    } else {
      sims[i] <- 1 / (1 + abs(trait_values[i] - trait_values[j]))
    }
  }
  
  # 5. Return mean similarity, handling NA values
  # If there are less than 2 valid similarities, return NA
  valid_sims <- sims[!is.na(sims)]
  if (length(valid_sims) < 2) {
    return(NA)
  }
  
  return(mean(valid_sims, na.rm = TRUE))
}

#' Build a model-specific covariance matrix for a subtree
#'
#' @param tree A subtree as a phylo object
#' @param model_name Model name (e.g. "Lambda1.0", "Lambda0.7", "OU_alpha1")
#' @return Covariance matrix induced by the tree and model
build_model_covariance_matrix <- function(tree, model_name) {
  n_tips <- length(tree$tip.label)
  if (n_tips < 2) {
    return(matrix(NA_real_, nrow = n_tips, ncol = n_tips,
                  dimnames = list(tree$tip.label, tree$tip.label)))
  }

  # Lambda models: reuse the same branch-length transformation used in simulation,
  # then compute the BM covariance matrix on the transformed subtree.
  if (grepl("^Lambda", model_name)) {
    lambda_value <- suppressWarnings(as.numeric(sub("^Lambda", "", model_name)))
    if (is.na(lambda_value) || lambda_value < 0) {
      stop(paste("Cannot parse a valid lambda value from model name:", model_name))
    }

    lambda_tree <- tree
    if (lambda_value < 1.0) {
      lambda_tree$edge.length <- lambda_tree$edge.length * lambda_value
      lambda_tree$edge.length <- lambda_tree$edge.length + (1 - lambda_value) * mean(lambda_tree$edge.length)
    }

    return(ape::vcv.phylo(lambda_tree))
  }

  # OU models: use the OU covariance function based on root-to-tip times,
  # pairwise patristic distances, and shared ancestry times.
  if (grepl("^OU_alpha", model_name)) {
    alpha_value <- suppressWarnings(as.numeric(sub("^OU_alpha", "", model_name)))
    if (is.na(alpha_value) || alpha_value <= 0) {
      stop(paste("Cannot parse a valid positive OU alpha from model name:", model_name))
    }

    node_depths <- ape::node.depth.edgelength(tree)
    tip_depths <- node_depths[seq_len(n_tips)]
    names(tip_depths) <- tree$tip.label

    mrca_matrix <- ape::mrca(tree)
    mrca_tips <- mrca_matrix[tree$tip.label, tree$tip.label, drop = FALSE]
    shared_times <- matrix(node_depths[mrca_tips],
                           nrow = n_tips,
                           ncol = n_tips,
                           dimnames = list(tree$tip.label, tree$tip.label))

    dist_matrix <- cophenetic.phylo(tree)
    cov_matrix <- (1 / (2 * alpha_value)) * exp(-alpha_value * dist_matrix) *
      (1 - exp(-2 * alpha_value * shared_times))
    diag(cov_matrix) <- (1 / (2 * alpha_value)) * (1 - exp(-2 * alpha_value * tip_depths))

    return(cov_matrix)
  }

  # Fallback: standard BM covariance.
  ape::vcv.phylo(tree)
}

#' Convert a covariance matrix to a correlation matrix
#'
#' @param cov_matrix Covariance matrix
#' @return Correlation matrix
cov_to_cor_matrix <- function(cov_matrix) {
  if (is.null(cov_matrix) || length(cov_matrix) == 0) {
    return(cov_matrix)
  }

  sd_vec <- sqrt(diag(cov_matrix))
  if (any(!is.finite(sd_vec)) || any(sd_vec <= 0)) {
    stop("Covariance matrix has non-positive or non-finite diagonal values")
  }

  cor_matrix <- cov_matrix / outer(sd_vec, sd_vec)
  diag(cor_matrix) <- 1
  return(cor_matrix)
}

#' Calculate tree-induced dependence metrics for a subtree
#'
#' @param tree A subtree as a phylo object
#' @param model_name Model name (e.g. "Lambda1.0", "OU_alpha1")
#' @return Named list with dependence diagnostics
calc_tree_dependence_metrics <- function(tree, model_name) {
  n_tips <- length(tree$tip.label)
  empty_result <- list(
    mean_offdiag_cor = NA_real_,
    max_offdiag_cor = NA_real_,
    ess_1 = NA_real_,   # mean-estimation ESS: 1' R^{-1} 1
    ess_2 = NA_real_    # spectral ESS
  )

  if (n_tips < 2) {
    return(empty_result)
  }

  cov_matrix <- build_model_covariance_matrix(tree, model_name)
  cor_matrix <- cov_to_cor_matrix(cov_matrix)

  # Force exact symmetry to reduce numerical noise
  cor_matrix <- (cor_matrix + t(cor_matrix)) / 2

  offdiag_vals <- cor_matrix[upper.tri(cor_matrix, diag = FALSE)]

  if (length(offdiag_vals) == 0 || any(!is.finite(offdiag_vals))) {
    return(empty_result)
  }

  # ---- ESS 1: mean-estimation ESS ----
  # n_eff,mean = 1' R^{-1} 1
  one_vec <- rep(1, n_tips)

  ess_1_value <- tryCatch({
    # Prefer Cholesky when possible
    chol_R <- chol(cor_matrix)
    inv_R <- chol2inv(chol_R)
    as.numeric(t(one_vec) %*% inv_R %*% one_vec)
  }, error = function(e) {
    # Fall back to a more general solver if Cholesky fails
    tryCatch({
      x <- qr.solve(cor_matrix, one_vec)
      as.numeric(crossprod(one_vec, x))
    }, error = function(e2) {
      NA_real_
    })
  })

  # ---- ESS 2: spectral ESS ----
  # n_eff,spec = n^2 / sum_{i,j} R_ij^2
  ess2_denom <- sum(cor_matrix * cor_matrix)
  ess_2_value <- if (is.finite(ess2_denom) && ess2_denom > 0) {
    (n_tips^2) / ess2_denom
  } else {
    NA_real_
  }

  list(
    mean_offdiag_cor = mean(offdiag_vals),
    max_offdiag_cor = max(offdiag_vals),
    ess_1 = ess_1_value,
    ess_2 = ess_2_value
  )
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

#' Calculate signal metrics for simulated traits with timing
#'
#' @param simulation_results Results from simulate_traits_for_subsets
#' @param tree Full tree
#' @param n_reps_to_analyze Number of replicates to analyze (for performance)
#' @param record_timing Whether to record computation time for each metric (default: TRUE)
#' @return A list containing:
#'   - signal_metrics: A comprehensive data frame with signal metrics
#'   - timing_stats: A data frame with timing statistics for each metric (if record_timing = TRUE)
analyze_simulated_traits_signal <- function(simulation_results, tree, n_reps_to_analyze = 100, record_timing = TRUE) {
  cat("Analyzing phylogenetic signal in simulated traits...\n")
  
  if (record_timing) {
    cat("  Timing recording is ENABLED\n")
  }

  all_results <- data.frame()
  
  # Initialize timing data structure if recording is enabled
  timing_data <- list()
  if (record_timing) {
    timing_data <- list(
      K = numeric(0),
      Lambda = numeric(0),
      MoransI = numeric(0),
      mpnns = numeric(0),
      dependence_metrics = numeric(0)
    )
  }

  # Get subset information from metadata
  subset_names <- simulation_results$metadata$subset_names
  subsets <- simulation_results$metadata$subsets

  if (is.null(subsets)) {
    warning("No subsets found in simulation results metadata. Cannot calculate signal metrics.")
    return(list(signal_metrics = all_results, timing_stats = NULL))
  }

  # Precompute phylogenetic distance matrix for the full tree
  # This will be used for mpnns calculation
  full_dist_matrix <- cophenetic.phylo(tree)

  process_model_results <- function(model_results, model_name) {
    model_df <- data.frame()

    for (subset_name in names(model_results)) {
      cat("  Processing", model_name, "for subset:", subset_name, "\n")

      subset_idx <- which(subset_names == subset_name)
      if (length(subset_idx) == 0) {
        warning(paste("Subset", subset_name, "not found in metadata"))
        next
      }

      subset_tips <- subsets[[subset_idx]]
      subset_traits <- model_results[[subset_name]]
      n_reps <- min(n_reps_to_analyze, ncol(subset_traits))
      subset_tree <- keep.tip(tree, subset_tips)

      # Time dependence metrics calculation
      dependence_metrics <- tryCatch({
        if (record_timing) {
          start_time <- proc.time()[3]
          result <- calc_tree_dependence_metrics(subset_tree, model_name)
          end_time <- proc.time()[3]
          timing_data$dependence_metrics <<- c(timing_data$dependence_metrics, end_time - start_time)
          result
        } else {
          calc_tree_dependence_metrics(subset_tree, model_name)
        }
      }, error = function(e) {
        cat("    Error calculating dependence metrics for", model_name, "subset", subset_name, ":", e$message, "\n")
        list(
          mean_offdiag_cor = NA_real_,
          max_offdiag_cor = NA_real_,
          ess_1 = NA_real_,
          ess_2 = NA_real_
        )
      })

      for (rep in 1:n_reps) {
        trait_values <- subset_traits[, rep]
        names(trait_values) <- subset_tips

        # Time K calculation
        K_value <- tryCatch({
          if (record_timing) {
            start_time <- proc.time()[3]
            result <- calc_blombergs_K(trait_values, subset_tree)
            end_time <- proc.time()[3]
            timing_data$K <<- c(timing_data$K, end_time - start_time)
            result
          } else {
            calc_blombergs_K(trait_values, subset_tree)
          }
        }, error = function(e) {
          cat("    Error calculating K for replicate", rep, ":", e$message, "\n")
          NA
        })

        # Time Lambda calculation
        lambda_value <- tryCatch({
          if (record_timing) {
            start_time <- proc.time()[3]
            result <- calc_pagels_lambda(trait_values, subset_tree)
            end_time <- proc.time()[3]
            timing_data$Lambda <<- c(timing_data$Lambda, end_time - start_time)
            result
          } else {
            calc_pagels_lambda(trait_values, subset_tree)
          }
        }, error = function(e) {
          cat("    Error calculating Lambda for replicate", rep, ":", e$message, "\n")
          NA
        })

        # Time Moran's I calculation
        moran_result <- tryCatch({
          if (record_timing) {
            start_time <- proc.time()[3]
            result <- calc_morans_I(trait_values, subset_tree)
            end_time <- proc.time()[3]
            timing_data$MoransI <<- c(timing_data$MoransI, end_time - start_time)
            result
          } else {
            calc_morans_I(trait_values, subset_tree)
          }
        }, error = function(e) {
          cat("    Error calculating Moran's I for replicate", rep, ":", e$message, "\n")
          list(observed = NA, expected = NA, sd = NA, p.value = NA)
        })

        # Time mpnns calculation
        mpnns_value <- tryCatch({
          if (record_timing) {
            start_time <- proc.time()[3]
            result <- calc_mpnns(subset_tips, trait_values, full_dist_matrix)
            end_time <- proc.time()[3]
            timing_data$mpnns <<- c(timing_data$mpnns, end_time - start_time)
            result
          } else {
            calc_mpnns(subset_tips, trait_values, full_dist_matrix)
          }
        }, error = function(e) {
          cat("    Error calculating mpnns for replicate", rep, ":", e$message, "\n")
          NA
        })

        model_df <- rbind(model_df, data.frame(
          Model = model_name,
          Subset = subset_name,
          Replicate = rep,
          K = K_value,
          Lambda = lambda_value,
          MoransI = moran_result$observed,
          MoransI_expected = moran_result$expected,
          MoransI_sd = moran_result$sd,
          MoransI_p = moran_result$p.value,
          mpnns = mpnns_value,
          mean_offdiag_cor = dependence_metrics$mean_offdiag_cor,
          max_offdiag_cor = dependence_metrics$max_offdiag_cor,
          ess_1 = dependence_metrics$ess_1,
          ess_2 = dependence_metrics$ess_2,
          stringsAsFactors = FALSE
        ))
      }
    }

    model_df
  }

  if (!is.null(simulation_results$Lambda)) {
    lambda_models <- names(simulation_results$Lambda)
    for (lambda_model in lambda_models) {
      all_results <- rbind(all_results, process_model_results(simulation_results$Lambda[[lambda_model]], lambda_model))
    }
  }

  if (!is.null(simulation_results$OU)) {
    ou_models <- names(simulation_results$OU)
    for (ou_model in ou_models) {
      all_results <- rbind(all_results, process_model_results(simulation_results$OU[[ou_model]], ou_model))
    }
  }

  # Calculate timing statistics if recording was enabled
  timing_stats <- NULL
  if (record_timing && length(timing_data) > 0) {
    cat("\nCalculating timing statistics...\n")
    
    timing_stats <- data.frame()
    
    for (metric_name in names(timing_data)) {
      times <- timing_data[[metric_name]]
      if (length(times) > 0) {
        timing_stats <- rbind(timing_stats, data.frame(
          Metric = metric_name,
          N = length(times),
          Total_Time = sum(times),
          Mean_Time = mean(times),
          SD_Time = sd(times),
          Min_Time = min(times),
          Max_Time = max(times),
          Median_Time = median(times),
          Q1_Time = quantile(times, 0.25),
          Q3_Time = quantile(times, 0.75),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Print summary
    cat("\nTiming Summary (seconds):\n")
    cat(strrep("=", 60), "\n")
    for (i in 1:nrow(timing_stats)) {
      row <- timing_stats[i, ]
      cat(sprintf("%-20s: N=%d, Total=%.3fs, Mean=%.3fs ± %.3fs, Min=%.3fs, Max=%.3fs\n",
                  row$Metric, row$N, row$Total_Time, row$Mean_Time, row$SD_Time,
                  row$Min_Time, row$Max_Time))
    }
    cat(strrep("=", 60), "\n")
    
    # Calculate overall statistics
    if (nrow(timing_stats) > 0) {
      total_calculations <- sum(timing_stats$N)
      total_time <- sum(timing_stats$Total_Time)
      cat(sprintf("\nOverall: %d calculations, Total time: %.3fs, Average per calculation: %.3fs\n",
                  total_calculations, total_time, total_time/total_calculations))
    }
  }

  return(list(
    signal_metrics = all_results,
    timing_stats = timing_stats
  ))
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
