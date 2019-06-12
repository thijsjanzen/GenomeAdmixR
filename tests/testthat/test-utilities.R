context("utilities")

test_that("utilities", {

  vx <- simulate_admixture(total_runtime = 5, seed = 1, number_of_founders = 50)
  testthat::expect_silent(
    plot_chromosome(vx$population[[1]]$chromosome1)
  )

  vx <- simulate_admixture(total_runtime = 100, seed = 1, number_of_founders = 2,
                           markers = seq(0,1,by=0.01))
  testthat::expect_silent(
    plot_over_time(vx$frequencies, focal_location = 0.5)
  )
})