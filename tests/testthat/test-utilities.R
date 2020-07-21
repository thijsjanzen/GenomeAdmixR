context("utilities")

test_that("utilities", {

  vx <- simulate_admixture(total_runtime = 5, seed = 1, number_of_founders = 50)
  testthat::expect_silent(
    plot_chromosome(vx$population[[1]]$chromosome1)
  )

  vx <- simulate_admixture(total_runtime = 100,
                           seed = 1,
                           number_of_founders = 2,
                           markers = seq(0, 1, by = 0.01))
  testthat::expect_silent(
    plot_over_time(vx$frequencies, focal_location = 0.5)
  )

  vy <- simulate_admixture_migration(migration_rate = 0.01, seed = 4)
  testthat::expect_error(plot_over_time(vy$frequencies, focal_location = 0.5))

  vy <- simulate_admixture_migration(migration_rate = 0.01, seed = 4,
                                     markers = 0.5)

  testthat::expect_silent(
    plot_over_time(vy$frequencies, focal_location = 0.5)
  )
})

test_that("initial_frequencies", {
  vx <- simulate_admixture_migration(total_runtime = 5, seed = 1,
                                     initial_frequencies = c(1, 1))

  vx <- simulate_admixture_migration(total_runtime = 5, seed = 1,
                           initial_frequencies = c(1, 1, 0, 0, 0, 0, 1, 1))


  vy <- simulate_admixture_migration(total_runtime = 5, seed = 1,
                                     initial_frequencies = list(c(1, 1, 0, 0),
                                                                c(0, 0, 1, 1)))

  testthat::expect_true(all.equal(vx, vy))

  testthat::expect_error(
    simulate_admixture_migration(total_runtime = 5, seed = 1,
                                 initial_frequencies = c(1, 1, 0, 0,
                                                         0, 0, 1, 1, 1))
  )


})