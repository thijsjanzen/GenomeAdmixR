context("simulate_Admixture")

test_that("simulate_admixture", {
  select_matrix <- matrix(NA, nrow = 2, ncol = 5)

  s <- 0.1
  select_matrix[1, ] <- c(0.5, 0.5, 0.5 + 0.5 * s, 0.5 + s, 0)
  select_matrix[2, ] <- c(0.6, 0.5, 0.5 + 0.5 * s, 0.5 + s, 0)

  testthat::expect_message(
    vx <- simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 1000,
                           morgan = 1,
                           select_matrix = select_matrix,
                           multiplicative_selection = FALSE)
  )
})


test_that("simulate admixture use", {
  testthat::expect_output(
    vx <- simulate_admixture(pop_size = 100,
                             number_of_founders = 2,
                             total_runtime = 100,
                             morgan = 1,
                             select_matrix = NA,
                             verbose = TRUE,
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

  testthat::expect_silent(
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


test_that("simulate admixture use, junctions", {
 # skip("simulate admixture use, junctions")
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

  testthat::expect_silent(
  GenomeAdmixR::plot_difference_frequencies(pop)
)


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

test_that("simulate admixture use, pop size", {
 # skip("simulate admixture use, pop size")
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
