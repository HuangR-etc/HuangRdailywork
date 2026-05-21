# 11_sanity_tree_helpers.R
# Helper functions for Figure 1 idealized 32-tip sanity check trees

#' Make an ultrametric balanced 32-tip tree
#'
#' @return Ultrametric phylo object with 32 tips labeled B001..B032
make_ultrametric_balanced_tree_32 <- function() {
  tree <- ape::stree(32, type = "balanced")
  tree$tip.label <- paste0("B", sprintf("%03d", seq_len(32)))
  tree$edge.length <- rep(1, nrow(tree$edge))
  tree$edge.length <- tree$edge.length /
    max(ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))])
  tree
}

#' Make an ultrametric ladder (pectinate) tree
#'
#' @param n_tips Number of tips (default 32)
#' @return Ultrametric phylo object
make_ultrametric_ladder_tree <- function(n_tips = 32) {
  merge_mat <- matrix(NA_integer_, nrow = n_tips - 1, ncol = 2)
  heights <- numeric(n_tips - 1)
  
  merge_mat[1, ] <- c(-1, -2)
  heights[1] <- 1
  
  if (n_tips >= 3) {
    for (i in 2:(n_tips - 1)) {
      merge_mat[i, ] <- c(i - 1, -(i + 1))
      heights[i] <- i
    }
  }
  
  hc <- list(
    merge = merge_mat,
    height = heights,
    order = seq_len(n_tips),
    labels = paste0("L", sprintf("%03d", seq_len(n_tips))),
    method = "complete",
    call = match.call()
  )
  class(hc) <- "hclust"
  
  tr <- ape::as.phylo(hc)
  tr$edge.length <- tr$edge.length /
    max(ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))])
  
  tr
}

#' Verify that a tree is 32-tip ultrametric
#'
#' @param tree Phylo object
#' @param tree_name Name for error messages
check_sanity_tree <- function(tree, tree_name = "tree") {
  stopifnot(ape::Ntip(tree) == 32)
  stopifnot(ape::is.ultrametric(tree, tol = 1e-8))
  invisible(TRUE)
}
