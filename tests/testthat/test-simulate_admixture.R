context("simulate_Admixture")

test_that("simulate admixture use, pop size", {
  pop <- simulate_admixture(pop_size = 100,
                            total_runtime = 3)

  testthat::expect_equal(length(pop$population), 100)

  in_pop <- list(pop$population[[1]],
                 pop$population[[2]])

  pop2 <- simulate_admixture(input_population = in_pop,
                             total_runtime = 3)

  testthat::expect_equal(length(pop2$population),
                         length(in_pop))

  pop3 <- simulate_admixture(input_population = in_pop,
                             pop_size = 100,
                             total_runtime = 3)

  testthat::expect_equal(length(pop3$population), 100)
})
