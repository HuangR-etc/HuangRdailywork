test_that("distance metrics and dispersed selection work", {
  tree <- ape::rtree(20)
  d <- patristic_matrix(tree)
  subset <- tree$tip.label[1:4]

  metrics <- distance_metrics(d, subset)
  expect_named(metrics, c("MinPD", "MeanPD", "MeanNND", "MaxPD"))
  expect_true(metrics$MaxPD >= metrics$MinPD)

  disp <- select_dispersed(dist_mat = d, candidates = tree$tip.label, size = 5)
  expect_length(disp$selected, 5)
  expect_true(disp$metrics$MinPD >= 0)
})

test_that("clustered main and fast methods both return valid subsets", {
  tree <- ape::rtree(24)
  d <- patristic_matrix(tree)

  main <- select_clustered(
    dist_mat = d,
    candidates = tree$tip.label,
    size = 6,
    method = "multistart_exchange"
  )
  fast <- select_clustered(
    dist_mat = d,
    candidates = tree$tip.label,
    size = 6,
    method = "fast"
  )

  expect_length(main$selected, 6)
  expect_length(fast$selected, 6)
  expect_true(main$metrics$MeanPD >= 0)
  expect_true(fast$metrics$MeanPD >= 0)
})
