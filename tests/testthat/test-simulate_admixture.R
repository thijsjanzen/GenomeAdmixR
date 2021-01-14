context("simulate_Admixture")



test_that("simulate admixture use, junctions", {
  markers <- 0.5
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           markers = markers,
                           track_junctions = TRUE)

  num_j <- length(vx$junctions)
  testthat::expect_gt(num_j, 0)
  testthat::expect_equal(num_j, 100)
})