context("simulate_Admixture")

test_that("simulate_admixture", {
  select_matrix <- matrix(NA, nrow = 2, ncol = 5)

  s <- 0.1
  select_matrix[1, ] <- c(0.5, 0.5, 0.5 + 0.5 * s, 0.5 + s, 0)
  select_matrix[2, ] <- c(0.6, 0.5, 0.5 + 0.5 * s, 0.5 + s, 0)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 1000,
                           morgan = 1,
                           seed = 42,
                           select_matrix = select_matrix,
                           multiplicative_selection = FALSE)

  vy <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 1000,
                           morgan = 1,
                           seed = 42,
                           select_matrix = select_matrix,
                           multiplicative_selection = TRUE)
  are_equal <- all.equal(vx$population, vy$population)
  if (length(are_equal) > 1) are_equal <- FALSE

  testthat::expect_true(!are_equal)
})


test_that("simulate admixture use", {

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = FALSE,
                           multiplicative_selection = TRUE)

  select_matrix <- matrix(NA, nrow = 1, ncol = 5)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            seed = 42,
                                            select_matrix = select_matrix,
                                            progress_bar = TRUE,
                                            track_junctions = FALSE,
                                            multiplicative_selection = TRUE))

  select_matrix <- matrix(NA, nrow = 1, ncol = 3)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            seed = 42,
                                            select_matrix = select_matrix,
                                            progress_bar = TRUE,
                                            track_junctions = FALSE,
                                            multiplicative_selection = TRUE))

  markers <- seq(from = 0.4, to = 0.6, length.out = 100)
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = FALSE,
                           markers = markers,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = TRUE,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = TRUE,
                           markers = markers,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           initial_frequencies = c(0.5, 0.5),
                           total_runtime = 100,
                           seed = 42)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           initial_frequencies = c(0.5, 0.6),
                           total_runtime = 100,
                           seed = 42)

  markers <- 0.5
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           seed = 42,
                           markers = markers)

  # code coverage for displaying functions:
  vx$population
  vx$population[[1]]
  plot(vx$population[[1]])
  print(vx$population)
  print(vx$population[[1]])
})


test_that("simulate admixture use, junctions", {
  markers <- 0.5
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           seed = 42,
                           markers = markers,
                           track_junctions = TRUE)
})