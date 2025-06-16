test_that("joyplot", {
  markers <- seq(from = 0.2, to = 0.3, length.out = 100)

  selected_pop <- simulate_admixture(module = ancestry_module(number_of_founders = 3,
                                                              morgan = 1,
                                                              markers = markers
                                                            ),
                                      pop_size = 1000,
                                      total_runtime = 101)

  testthat::expect_lte(max(selected_pop$frequencies$location), max(markers))
  testthat::expect_gte(min(selected_pop$frequencies$location), min(markers))


  px <- plot_joyplot_frequencies(selected_pop$frequencies,
                                 time_points = c(0, 50, 100), picked_ancestor = "ALL")
  testthat::expect_identical(px$labels$x, "Location (Morgan)")
  testthat::expect_identical(px$labels$y, "Time")
  testthat::expect_identical(px$labels$fill, "Ancestor")

  px <- plot_joyplot_frequencies(selected_pop$frequencies,
                                 time_points = c(0, 50, 100), picked_ancestor = 1)
  testthat::expect_identical(px$labels$x, "Location (Morgan)")
  testthat::expect_identical(px$labels$y, "Time")
  testthat::expect_identical(px$labels$fill, "Ancestor")
  ax <- as.numeric(unique(px$data$ancestor)[[1]])
  testthat::expect_identical(ax, 2)
  # 2 because 1 is the second in levels [0, 1, 2]
})
