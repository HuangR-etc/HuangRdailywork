context("Distance metrics")

test_that("distance_metrics returns correct structure", {
  library(ape)
  tree <- rtree(10)
  D <- ape::cophenetic.phylo(tree)
  subset <- tree$tip.label[1:4]

  m <- distance_metrics(D, subset)
  expect_type(m, "list")
  expect_named(m, c("MinPD", "MaxPD", "MeanPD", "MeanNND"))
  expect_true(m$MinPD >= 0)
  expect_true(m$MaxPD >= m$MinPD)
})

test_that("single species returns zeros", {
  D <- matrix(c(0), nrow = 1, dimnames = list("sp1", "sp1"))
  m <- distance_metrics(D, "sp1")
  expect_equal(m$MinPD, 0)
  expect_equal(m$MaxPD, 0)
  expect_equal(m$MeanPD, 0)
  expect_equal(m$MeanNND, 0)
})

test_that("individual metric functions work", {
  library(ape)
  tree <- rtree(10)
  D <- ape::cophenetic.phylo(tree)
  subset <- tree$tip.label[1:4]

  expect_true(min_pd(D, subset) >= 0)
  expect_true(max_pd(D, subset) >= min_pd(D, subset))
  expect_true(mean_pd(D, subset) >= 0)
  expect_true(mean_nnd(D, subset) >= 0)
})

test_that("nearest_neighbor_distances returns named vector", {
  library(ape)
  tree <- rtree(10)
  D <- ape::cophenetic.phylo(tree)
  subset <- tree$tip.label[1:4]

  nnd <- nearest_neighbor_distances(D, subset)
  expect_type(nnd, "double")
  expect_named(nnd, subset)
})
