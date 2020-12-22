context("tajima")

testthat::test_that("tajima", {
  pop <- simulate_admixture(pop_size = 100,
                            number_of_founders = 2,
                            seed = 666,
                            total_runtime = 10)

  t <- 10

  found <- calculate_tajima_d(pop$population)$D

  while (t < 1000) {
    pop <- simulate_admixture(pop$population,
                              pop_size = 100,
                              seed = t,
                              total_runtime = 10,
                              progress_bar = FALSE)
    found <- c(found,
               calculate_tajima_d(pop$population)$D
    )
    t <- t + 10
  }

  testthat::expect_true(mean(found, na.rm = T) < 2)
})

testthat::test_that("tajima abuse", {
  pop <- simulate_admixture(pop_size = 100,
                            number_of_founders = 2,
                            seed = 666,
                            total_runtime = 10)

  testthat::expect_error(calculate_tajima_d(pop$initial_frequency))
})
