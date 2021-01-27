context("test input data")

test_that("input data", {
  testthat::skip_on_os("solaris")
  message("test input data")

  chosen_markers <- 1:100

  fake_input_data1 <- list()
  fake_input_data1$genomes <- matrix(data = 1, nrow = 100, ncol = 100)
  fake_input_data1$markers <- chosen_markers

  fake_input_data2 <- list()
  fake_input_data2$genomes <- matrix(data = 2, nrow = 100, ncol = 100)
  fake_input_data2$markers <- chosen_markers

  combined_data <- combine_input_data(input_data_list = list(fake_input_data1,
                                                             fake_input_data2),
                                      frequencies = c(0.3, 0.7),
                                      pop_size = 1000)

  vv <- table(combined_data$genomes)
  a <- vv[[1]] / (sum(vv))
  testthat::expect_equal(a, 0.3, tolerance = 0.1)

  testthat::expect_error(
    combined_data <- combine_input_data(input_data_list = list(fake_input_data1,
                                                               fake_input_data2),
                                        frequencies = c(0.1, 0.7, 0.2),
                                        pop_size = 1000)
  )

})