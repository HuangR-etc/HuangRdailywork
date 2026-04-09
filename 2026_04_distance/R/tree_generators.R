# Tree generation functions for phylogenetic dispersed subset analysis
# This module generates balanced and ladder (unbalanced) trees

library(ape)
library(phytools)

#' Generate a balanced tree with given number of tips
#'
#' @param n_tips Number of tips in the tree
#' @param branch_length_mean Mean branch length (default: 1.0)
#' @param branch_length_sd Standard deviation of branch length (default: 0.2)
#' @param tip_prefix Prefix for tip labels (default: "sp")
#' @return A phylo object
generate_balanced_tree <- function(n_tips, branch_length_mean = 1.0, branch_length_sd = 0.2, tip_prefix = "sp") {
  # Create a tree with random topology and branch lengths
  tree <- rtree(n_tips, br = branch_length_mean, rooted = TRUE)
  
  # Make it more balanced by reordering
  # We'll use a balanced tree structure from stree and transfer branch lengths
  if (n_tips > 1) {
    # Create a balanced tree structure
    balanced_structure <- stree(n_tips, type = "balanced")
    
    # Keep the branch lengths from the random tree
    # We need to ensure both trees have the same number of edges
    if (length(tree$edge.length) == length(balanced_structure$edge.length)) {
      balanced_structure$edge.length <- tree$edge.length
      tree <- balanced_structure
    }
    # If edge counts don't match, keep the random tree (it will still work)
  }
  
  # Assign tip labels
  tree$tip.label <- paste0(tip_prefix, sprintf(paste0("%0", nchar(as.character(n_tips)), "d"), 1:n_tips))
  
  return(tree)
}

#' Generate a ladder (unbalanced) tree with given number of tips
#'
#' @param n_tips Number of tips in the tree
#' @param branch_length_mean Mean branch length (default: 1.0)
#' @param branch_length_sd Standard deviation of branch length (default: 0.2)
#' @param tip_prefix Prefix for tip labels (default: "sp")
#' @return A phylo object
generate_ladder_tree <- function(n_tips, branch_length_mean = 1.0, branch_length_sd = 0.2, tip_prefix = "sp") {
  # Create an unbalanced (ladder) tree
  # We'll create a pectinate tree (completely unbalanced)
  tree <- stree(n_tips, type = "left")
  
  # Generate random branch lengths
  n_edges <- nrow(tree$edge)
  tree$edge.length <- abs(rnorm(n_edges, mean = branch_length_mean, sd = branch_length_sd))
  
  # Ensure no zero or negative branch lengths
  tree$edge.length[tree$edge.length <= 0] <- 0.1
  
  # Assign tip labels
  tree$tip.label <- paste0(tip_prefix, sprintf(paste0("%0", nchar(as.character(n_tips)), "d"), 1:n_tips))
  
  return(tree)
}

#' Generate all trees needed for the analysis
#'
#' @param cfg Configuration list
#' @return A list containing all trees
generate_all_trees <- function(cfg) {
  cat("Generating trees...\n")
  
  trees <- list()
  
  # Generate large balanced tree
  cat("  Generating large balanced tree (n =", cfg$large_n, ")...\n")
  trees$large_balanced <- generate_balanced_tree(cfg$large_n, tip_prefix = "LB")
  
  # Generate large ladder (unbalanced) tree
  cat("  Generating large ladder tree (n =", cfg$large_n, ")...\n")
  trees$large_ladder <- generate_ladder_tree(cfg$large_n, tip_prefix = "LL")
  
  # Generate small balanced tree
  cat("  Generating small balanced tree (n =", cfg$small_n, ")...\n")
  trees$small_balanced <- generate_balanced_tree(cfg$small_n, tip_prefix = "SB")
  
  # If tree_reps > 1, generate replicates
  if (cfg$tree_reps > 1) {
    cat("  Generating", cfg$tree_reps, "replicates for each tree type...\n")
    
    # Create lists for replicates
    trees$large_balanced_reps <- list()
    trees$large_ladder_reps <- list()
    trees$small_balanced_reps <- list()
    
    for (i in 1:cfg$tree_reps) {
      trees$large_balanced_reps[[i]] <- generate_balanced_tree(
        cfg$large_n, 
        tip_prefix = paste0("LB", i, "_")
      )
      trees$large_ladder_reps[[i]] <- generate_ladder_tree(
        cfg$large_n, 
        tip_prefix = paste0("LL", i, "_")
      )
      trees$small_balanced_reps[[i]] <- generate_balanced_tree(
        cfg$small_n, 
        tip_prefix = paste0("SB", i, "_")
      )
    }
  }
  
  cat("Tree generation complete.\n")
  return(trees)
}

#' Plot and save tree visualizations
#'
#' @param trees List of trees
#' @param cfg Configuration list
plot_trees <- function(trees, cfg) {
  cat("Plotting trees...\n")
  
  # Create figure directory if it doesn't exist
  fig_dir <- cfg$figures_dir
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Plot large balanced tree
  pdf(file.path(fig_dir, "tree_large_balanced.pdf"), width = 10, height = 8)
  plot(trees$large_balanced, main = "Large Balanced Tree (n=256)", cex = 0.6)
  dev.off()
  
  # Plot large ladder tree
  pdf(file.path(fig_dir, "tree_large_ladder.pdf"), width = 10, height = 8)
  plot(trees$large_ladder, main = "Large Ladder Tree (n=256)", cex = 0.6)
  dev.off()
  
  # Plot small balanced tree
  pdf(file.path(fig_dir, "tree_small_balanced.pdf"), width = 8, height = 6)
  plot(trees$small_balanced, main = "Small Balanced Tree (n=32)", cex = 0.8)
  dev.off()
  
  cat("Tree plots saved to", fig_dir, "\n")
}

# Test function
if (sys.nframe() == 0) {
  # Load configuration
  source("../config/analysis_config.R")
  
  # Generate trees
  test_trees <- generate_all_trees(cfg)
  
  # Print tree summaries
  cat("\nTree summaries:\n")
  cat("  Large balanced tree:", length(test_trees$large_balanced$tip.label), "tips\n")
  cat("  Large ladder tree:", length(test_trees$large_ladder$tip.label), "tips\n")
  cat("  Small balanced tree:", length(test_trees$small_balanced$tip.label), "tips\n")
  
  # Plot trees
  plot_trees(test_trees, cfg)
}
