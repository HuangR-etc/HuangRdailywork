test_that("covariance helpers and MIESS behave on simple cases", {
  r <- diag(5)
  dep <- dependence_metrics(r)
  expect_equal(dep$MeanOffCor, 0)
  expect_equal(dep$MaxOffCor, 0)
  expect_equal(dep$MIESS, 5)

  tree <- ape::rtree(12)
  v_bm <- phylo_covariance(tree, model = "BM")
  v_lambda0 <- phylo_covariance(tree, model = "lambda", lambda = 0)
  v_eb0 <- phylo_covariance(tree, model = "EB", eb_rate = 0)
  v_eb_branch0 <- phylo_covariance(
    tree,
    model = "EB",
    eb_rate = 0,
    eb_method = "branch_transform"
  )
  expect_equal(dim(v_bm), c(12, 12))
  expect_true(all(v_lambda0[upper.tri(v_lambda0)] == 0))
  expect_equal(v_eb0, v_bm)
  expect_equal(v_eb_branch0, v_bm)
})

test_that("EB package default matches manuscript simple covariance", {
  tree <- ape::rtree(10)
  v_pkg <- phylo_covariance(tree, model = "EB", eb_rate = -0.5)
  v_simple <- make_eb_covariance_simple(tree, r = -0.5)
  v_branch <- phylo_covariance(
    tree,
    model = "EB",
    eb_rate = -0.5,
    eb_method = "branch_transform"
  )

  expect_equal(v_pkg, v_simple)
  expect_false(isTRUE(all.equal(v_pkg, v_branch)))
})

test_that("PIESS utilities return expected columns", {
  bench <- run_independent_benchmark_curve(4:6, n_sim = 200, seed = 10)
  expect_true(all(c("Metric", "Benchmark_N", "Interval_Width_95") %in% names(bench)))

  target <- bench[bench$Benchmark_N == 6 & bench$Metric == "RMSE", , drop = FALSE]
  ess <- estimate_prediction_metric_ess(target, bench)
  expect_true("Prediction_Metric_ESS_Label" %in% names(ess))
})

test_that("phylo_test_ess can compute PIESS under manuscript-default EB", {
  tree <- ape::rtree(24)
  res <- phylo_test_ess(
    tree = tree,
    candidates = tree$tip.label,
    size = 6,
    subset_type = "dispersed",
    n_random = 20,
    cov_model = "EB",
    eb_rate = -0.1,
    compute_piess = TRUE,
    piess_n_sim = 200,
    seed = 123
  )

  expect_s3_class(res, "phylo_test_ess_result")
  expect_true(all(c("MeanOffCor", "MaxOffCor", "MIESS") %in% names(res$dependence_metrics)))
  expect_true("Prediction_Metric_ESS_Label" %in% names(res$piess$piess))
})
