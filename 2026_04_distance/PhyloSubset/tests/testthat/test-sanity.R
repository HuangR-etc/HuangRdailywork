context("Sanity checks on balanced and ladder trees")

test_that("dispersed subset has higher MinPD than random on balanced tree", {
  library(ape)
  tree <- ape::stree(16, type = "balanced")
  tree$edge.length <- rep(1, ape::Nedge(tree))
  candidates <- tree$tip.label

  D <- patristic_matrix(tree, candidates)

  disp <- select_dispersed(dist_mat = D, candidates = candidates, size = 4,
                           refine = TRUE, seed = 42)

  rand_metrics <- random_baseline(D, candidates, size = 4, n = 100, seed = 99)

  p_val <- empirical_p_value(disp$metrics$MinPD, rand_metrics$MinPD, "greater")
  expect_true(p_val <= 0.05)
})

test_that("clustered subset has lower MeanPD than random on ladder tree", {
  library(ape)
  tree <- ape::stree(16, type = "left")
  tree$edge.length <- rep(1, ape::Nedge(tree))
  candidates <- tree$tip.label

  D <- patristic_matrix(tree, candidates)

  clust <- select_clustered(dist_mat = D, candidates = candidates, size = 4)

  rand_metrics <- random_baseline(D, candidates, size = 4, n = 100, seed = 99)

  p_val <- empirical_p_value(clust$metrics$MeanPD, rand_metrics$MeanPD, "less")
  expect_true(p_val <= 0.05)
})

test_that("input checks catch invalid tree", {
  expect_error(check_phylo_input(list()))
  expect_error(check_phylo_input("not_a_tree"))
})

test_that("input checks catch missing branch lengths", {
  library(ape)
  tree <- ape::stree(8, type = "balanced")
  expect_error(check_phylo_input(tree))
})

test_that("match_species handles underscore conversion", {
  library(ape)
  tree <- rtree(5)
  tree$tip.label <- paste0("Species_", 1:5)
  result <- match_species(tree, c("Species 1", "Species 6"), fuzzy = TRUE)
  expect_true("Species_1" %in% result$matched)
  expect_true("Species 6" %in% result$unmatched)
})
