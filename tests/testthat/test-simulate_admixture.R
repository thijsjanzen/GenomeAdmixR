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

