test_that("phylo_test_ess returns a structured result", {
  tree <- ape::rtree(30)
  res <- phylo_test_ess(
    tree = tree,
    candidates = tree$tip.label,
    size = 6,
    subset_type = "clustered",
    clustered_method = "fast",
    n_random = 20,
    seed = 123
  )

  expect_s3_class(res, "phylo_test_ess_result")
  expect_length(res$selected, 6)
  expect_true(all(c("MeanOffCor", "MaxOffCor", "MIESS") %in% names(res$dependence_metrics)))
  expect_equal(nrow(as.data.frame(res)), 1)
})
