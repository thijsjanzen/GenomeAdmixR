context("simulate_Admixture")

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

  select_matrix <- matrix(NA, nrow=1, ncol=5)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            seed = 42,
                                            select_matrix = select_matrix,
                                            progress_bar = TRUE,
                                            track_junctions = FALSE,
                                            multiplicative_selection = TRUE))

  select_matrix <- matrix(NA, nrow=1,ncol=3)
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


  vx <- simulate_admixture(pop_size = 100,
                                   number_of_founders = 2,
                                   total_runtime = 100,
                                   morgan = 1,
                                   seed = 42,
                                   select_matrix = NA,
                                   progress_bar = TRUE,
                                   track_junctions = TRUE,
                                   markers = markers,
                                   multiplicative_selection = FALSE)

  # code coverage for displaying functions:
  vx$population
  vx$population[[1]]
  plot(vx$population[[1]])
  print(vx$population)
  print(vx$population[[1]])

})