context("Random baseline and p-values")

test_that("random_subsets returns correct number of subsets", {
  candidates <- paste0("sp", 1:20)
  subsets <- random_subsets(candidates, size = 5, n = 10)
  expect_length(subsets, 10)
  expect_true(all(lengths(subsets) == 5))
})

test_that("random_subsets respects seed", {
  candidates <- paste0("sp", 1:20)
  s1 <- random_subsets(candidates, size = 5, n = 3, seed = 42)
  s2 <- random_subsets(candidates, size = 5, n = 3, seed = 42)
  expect_identical(s1, s2)
})

test_that("empirical_p_value returns values between 0 and 1", {
  p <- empirical_p_value(0.8, rnorm(100, 0.5, 0.2), "greater")
  expect_true(p >= 0 && p <= 1)
})

test_that("empirical_p_value with extreme values", {
  # Observed far above all random values
  p <- empirical_p_value(10, rnorm(100, 0, 1), "greater")
  expect_equal(p, 1 / 101, tolerance = 1e-6)

  # Observed far below all random values
  p <- empirical_p_value(-10, rnorm(100, 0, 1), "less")
  expect_equal(p, 1 / 101, tolerance = 1e-6)
})

test_that("compare_to_random returns correct structure", {
  library(ape)
  tree <- rtree(20)
  D <- ape::cophenetic.phylo(tree)
  candidates <- tree$tip.label
  obs <- distance_metrics(D, sample(candidates, 5))
  rand <- random_baseline(D, candidates, size = 5, n = 50)

  comp <- compare_to_random(obs, rand, type = "dispersed")
  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 4)
  expect_true(all(c("Metric", "Observed", "P_value", "SES") %in% names(comp)))
})
