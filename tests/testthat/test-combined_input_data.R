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

  chosen_markers <- 1:1000

  fake_input_data1 <- list()
  fake_input_data1$genomes <- matrix(data = 1, nrow = 100, ncol = 1000)
  fake_input_data1$markers <- chosen_markers

  fake_input_data2 <- list()
  fake_input_data2$genomes <- matrix(data = 2, nrow = 100, ncol = 1000)
  fake_input_data2$markers <- chosen_markers

  combined_data <- combine_input_data(input_data_list = list(fake_input_data1,
                                                             fake_input_data2),
                                      frequencies = c(0.5, 0.5),
                                      pop_size = 100)

  simulation_result <- simulate_admixture_data(input_data = combined_data,
                                               pop_size = 1000,
                                               total_runtime = 10,
                                               morgan = 1)

  num_j <- c()
  for (i in 1:length(simulation_result$population)) {
    focal_indiv <- simulation_result$population[[i]]
    chrom1 <- focal_indiv$chromosome1
    chrom2 <- focal_indiv$chromosome2
    num_j_1 <- sum(abs(diff(chrom1[,2])))
    num_j_2 <- sum(abs(diff(chrom2[,2])))
    num_j <- c(num_j, num_j_1, num_j_2)
  }
  expected_j <- junctions::number_of_junctions(N = 1000, R = 1000, t = 10)
  testthat::expect_equal(mean(num_j), expected_j, tolerance = 1)
})