context("Covariance and dependence metrics")

test_that("cov_to_cor returns valid correlation matrix", {
  V <- matrix(c(4, 1, 1, 9), nrow = 2)
  R <- cov_to_cor(V)
  expect_equal(diag(R), c(1, 1))
  expect_true(abs(R[1, 2]) <= 1)
})

test_that("dependence_metrics returns correct structure", {
  R <- diag(5)
  dep <- dependence_metrics(R)
  expect_type(dep, "list")
  expect_named(dep, c("MeanOffCor", "MaxOffCor", "MeanESS"))
  expect_equal(dep$MeanOffCor, 0)
  expect_equal(dep$MaxOffCor, 0)
  expect_equal(dep$MeanESS, 5)
})

test_that("mean_ess returns nominal size for identity matrix", {
  R <- diag(10)
  ess <- mean_ess(R)
  expect_equal(ess, 10, tolerance = 1e-8)
})

test_that("mean_ess handles singular matrix gracefully", {
  R <- matrix(1, nrow = 3, ncol = 3)
  ess <- mean_ess(R, method = "qr")
  expect_true(is.na(ess) || ess > 0)
})

test_that("phylo_covariance BM returns valid matrix", {
  library(ape)
  tree <- rtree(10)
  V <- phylo_covariance(tree, model = "BM")
  expect_equal(dim(V), c(10, 10))
  expect_true(all(diag(V) > 0))
})

test_that("lambda = 0 gives diagonal covariance", {
  library(ape)
  tree <- rtree(10)
  V <- phylo_covariance(tree, model = "lambda", lambda = 0)
  expect_true(all(V[upper.tri(V)] == 0))
})
