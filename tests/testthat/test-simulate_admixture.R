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
                           select_matrix = select_matrix,
                           multiplicative_selection = FALSE)
})


test_that("simulate admixture use", {

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = FALSE,
                           multiplicative_selection = TRUE)

  select_matrix <- matrix(NA, nrow = 1, ncol = 5)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            select_matrix = select_matrix,
                                            progress_bar = TRUE,
                                            track_junctions = FALSE,
                                            multiplicative_selection = TRUE))

  select_matrix <- matrix(NA, nrow = 1, ncol = 3)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            select_matrix = select_matrix,
                                            progress_bar = TRUE,
                                            track_junctions = FALSE,
                                            multiplicative_selection = TRUE))

  markers <- seq(from = 0.4, to = 0.6, length.out = 100)
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = FALSE,
                           markers = markers,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = TRUE,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = TRUE,
                           markers = markers,
                           multiplicative_selection = TRUE)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           initial_frequencies = c(0.5, 0.5),
                           total_runtime = 100)

  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           initial_frequencies = c(0.5, 0.6),
                           total_runtime = 100)

  markers <- 0.5
  vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
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
                           markers = markers,
                           track_junctions = TRUE)

  num_j <- length(vx$junctions)
  testthat::expect_gt(num_j, 0)
  testthat::expect_equal(num_j, 100)
})

test_that("simulate admixture use, markers", {
  pop <- simulate_admixture(pop_size = 1000,
                            number_of_founders = 2,
                            total_runtime = 3,
                            markers = seq(0, 1, length.out = 1000),
                            morgan = 1)

  a <- subset(pop$final_frequency, pop$final_frequency$location < 0.5)
  b <-  subset(pop$final_frequency, pop$final_frequency$location > 0.5)
  testthat::expect_equal(dim(a), dim(b))

  GenomeAdmixR::plot_difference_frequencies(pop)



  pop <- simulate_admixture(pop_size = 1000,
                            number_of_founders = 2,
                            total_runtime = 3,
                            markers = seq(0, 3, length.out = 1000),
                            morgan = 3)

  a <- subset(pop$final_frequency, pop$final_frequency$location < 1.5)
  b <-  subset(pop$final_frequency, pop$final_frequency$location > 1.5)
  testthat::expect_equal(dim(a), dim(b))

  GenomeAdmixR::plot_difference_frequencies(pop)

  pop <- simulate_admixture(pop_size = 1000,
                            number_of_founders = 2,
                            total_runtime = 3,
                            markers = seq(0, 5, length.out = 1000),
                            morgan = 5)

  a <- subset(pop$final_frequency, pop$final_frequency$location < 2.5)
  b <-  subset(pop$final_frequency, pop$final_frequency$location > 2.5)
  testthat::expect_equal(dim(a), dim(b))
  GenomeAdmixR::plot_difference_frequencies(pop)
})

test_that("simulate admixture use, threads", {
  skip_on_cran()
  population_size <- 1000

  t1 <- Sys.time()
  vx <- simulate_admixture(pop_size = population_size,
                           total_runtime = 1000,
                           morgan = 1,
                           num_threads = 1)

  t2 <- Sys.time()
  vy <- simulate_admixture(pop_size = population_size,
                           total_runtime = 1000,
                           morgan = 1,
                           num_threads = -1)
  t3 <- Sys.time()

  time_one_thread = difftime(t2, t1, units = "secs")[[1]]
  time_all_threads = difftime(t3, t2, units = "secs")[[1]]
  testthat::expect_lt(time_all_threads, time_one_thread)
})
