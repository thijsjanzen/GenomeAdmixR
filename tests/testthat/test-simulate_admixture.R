context("simulate_Admixture")


test_that("simulate admixture use", {
  testthat::expect_output(
    vx <- simulate_admixture(pop_size = 100,
                             number_of_founders = 2,
                             total_runtime = 100,
                             morgan = 1,
                             select_matrix = NA,
                             progress_bar = TRUE,
                             track_junctions = FALSE,
                             multiplicative_selection = TRUE)
  )

  select_matrix <- matrix(NA, nrow = 1, ncol = 5)
  testthat::expect_message(
    testthat::expect_error(simulate_admixture(pop_size = 100,
                                              number_of_founders = 2,
                                              total_runtime = 100,
                                              morgan = 1,
                                              select_matrix = select_matrix,
                                              track_junctions = FALSE,
                                              multiplicative_selection = TRUE))
  )
  select_matrix <- matrix(NA, nrow = 1, ncol = 3)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            select_matrix = select_matrix,
                                            track_junctions = FALSE,
                                            multiplicative_selection = TRUE))

  markers <- seq(from = 0.4, to = 0.6, length.out = 100)
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           select_matrix = NA,
                           track_junctions = FALSE,
                           markers = markers,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           select_matrix = NA,
                           track_junctions = TRUE,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           select_matrix = NA,
                           track_junctions = TRUE,
                           markers = markers,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           initial_frequencies = c(0.5, 0.5),
                           total_runtime = 100)

  testthat::expect_message(
    vx <- simulate_admixture(pop_size = 100,
                             number_of_founders = 2,
                             initial_frequencies = c(0.5, 0.6),
                             total_runtime = 100)
  )
  markers <- 0.5
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           markers = markers)

  # code coverage for displaying functions:
  testthat::expect_silent(plot(vx$population[[1]]))

  testthat::expect_output(
    testthat::expect_equal( print(vx$population),
                            "Population with 100 individuals")
  )
  testthat::expect_output(print(vx$population[[1]]))
})
