test_that("joyplot", {

  s <- 0.5
  select_matrix <- matrix(nrow = 1, ncol = 5)
  select_matrix[1, ] <- c(0.25, 1.0, 1 + 0.5 * s, 1 + s, 0)

  markers <- seq(from = 0.2, to = 0.3, length.out = 100)

  selected_pop <- simulate_admixture(pop_size = 1000,
                                     number_of_founders = 3,
                                     total_runtime = 100,
                                     morgan = 1,
                                     select_matrix,
                                     markers = markers)

  px <- plot_joyplot_frequencies(selected_pop$frequencies,
                           time_points = c(0, 50, 100), picked_ancestor = "ALL")
  testthat::expect_identical(px$labels$x, "Location (Morgan)")
  testthat::expect_identical(px$labels$y, "Time")
  testthat::expect_identical(px$labels$fill, "Ancestor")
  testthat::expect_identical(px$labels$height, "frequency")
})
