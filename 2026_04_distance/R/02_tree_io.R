# 02_tree_io.R
# Tree reading, cleaning, and taxonomic extraction

#' Read and clean the mammal tree
#'
#' Reads the full mammal tree from NEXUS file and drops specified anomalous tips.
#'
#' @param tree_file Path to the NEXUS tree file
#' @param drop_tips Character vector of tip labels to drop
#' @return Cleaned phylo object
read_clean_mammal_tree <- function(tree_file = TREE_FILE,
                                   drop_tips = c(
                                     "_Anolis_carolinensis",
                                     "Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA"
                                   )) {
  tree <- ape::read.nexus(tree_file)
  tree <- ape::drop.tip(tree, intersect(drop_tips, tree$tip.label))
  tree
}

#' Extract a taxonomic pool subtree
#'
#' Filters the taxonomy table to a given rank and value, then drops all tips
#' from the tree that are not in that pool.
#'
#' @param tree Full phylo object
#' @param tax Taxonomy data frame
#' @param rank_col Column name in tax for the rank (e.g. "fam", "order")
#' @param rank_value Value to filter by (e.g. "CRICETIDAE", "Carnivora")
#' @param require_extant If TRUE, filter to extant species only
#' @param require_sampled If TRUE, filter to sampled species only
#' @return Subtree phylo object
extract_taxonomic_pool_tree <- function(tree,
                                        tax,
                                        rank_col,
                                        rank_value,
                                        require_extant = TRUE,
                                        require_sampled = TRUE) {
  tax_sub <- tax[toupper(tax[[rank_col]]) == toupper(rank_value), ]
  
  if (require_extant && "extinct." %in% names(tax_sub)) {
    tax_sub <- tax_sub[tax_sub[["extinct."]] == 0, ]
  }
  
  if (require_sampled && "samp" %in% names(tax_sub)) {
    tax_sub <- tax_sub[tax_sub$samp == "sampled", ]
  }
  
  tips <- intersect(tax_sub$tiplabel, tree$tip.label)
  
  ape::drop.tip(tree, setdiff(tree$tip.label, tips))
}

#' Build nested candidate pools for sensitivity analysis
#'
#' From a pool tree, randomly sample the largest pool, then subsample
#' progressively smaller pools from it. Returns both tip name vectors
#' and corresponding subtrees.
#'
#' @param pool_tree The full pool tree (e.g. Cricetidae)
#' @param pool_sizes Vector of pool sizes (ascending)
#' @param seed Random seed for reproducibility
#' @return List with $pools (named list of tip vectors) and $pool_trees (named list of phylo)
build_nested_candidate_pools <- function(pool_tree,
                                         pool_sizes = c(32, 64, 128, 256, 512),
                                         seed = GLOBAL_SEED) {
  max_N <- max(pool_sizes)
  
  if (ape::Ntip(pool_tree) < max_N) {
    stop("Not enough tips in pool_tree for max pool size.")
  }
  
  set.seed(seed)
  
  pool_sizes_desc <- sort(pool_sizes, decreasing = TRUE)
  pools <- list()
  
  pools[[paste0("C", max_N)]] <- sample(pool_tree$tip.label, max_N)
  
  for (i in seq_len(length(pool_sizes_desc) - 1)) {
    parent_N <- pool_sizes_desc[i]
    child_N <- pool_sizes_desc[i + 1]
    
    parent_label <- paste0("C", parent_N)
    child_label <- paste0("C", child_N)
    
    pools[[child_label]] <- sample(pools[[parent_label]], child_N)
  }
  
  pool_sizes_asc <- sort(pool_sizes)
  pools <- pools[paste0("C", pool_sizes_asc)]
  
  pool_trees <- lapply(pools, function(tips) {
    ape::drop.tip(pool_tree, setdiff(pool_tree$tip.label, tips))
  })
  
  list(
    pools = pools,
    pool_trees = pool_trees
  )
}
