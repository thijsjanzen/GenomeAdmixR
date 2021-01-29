context("test simulate admixture data")

test_that("simulate_admixture_data", {
  testthat::skip_on_os("solaris")
  message("test simulate_admixture_data")

  num_markers <- 100
  num_indiv <- 100
  chosen_markers <- 1:num_markers

  fake_input_data1 <- list()
  fake_input_data1$genomes <- matrix(data = sample(x = 1:2,
                                                   size = num_indiv *
                                                           num_markers,
                                                   replace = T),
                                     nrow = num_indiv,
                                     ncol = num_markers)


  fake_input_data1$markers <- chosen_markers

  fake_input_data2 <- list()
  fake_input_data2$genomes <- matrix(data = sample(x = 3:4,
                                                   size = num_indiv *
                                                           num_markers,
                                                   replace = T),
                                     nrow = num_indiv,
                                     ncol = num_markers)
  fake_input_data2$markers <- chosen_markers

  combined_data <- combine_input_data(input_data_list = list(fake_input_data1,
                                                             fake_input_data2),
                                      frequencies = c(0.5, 0.5),
                                      pop_size = 1000)


  simul_pop <- simulate_admixture_data(input_data = combined_data,
                                       pop_size = 100,
                                       total_runtime = 100,
                                       markers = chosen_markers,
                                       morgan = 1,
                                       verbose = FALSE)

  testthat::expect_silent(
    plot_chromosome(simul_pop$population[[1]]$chromosome1)
  )

  testthat::expect_silent(
    plot_difference_frequencies(results = simul_pop)
  )
  testthat::expect_silent(
    calculate_allele_frequencies(source_pop = simul_pop,
                                 progress_bar = FALSE)

  )
  testthat::expect_silent(
    plot_frequencies(simul_pop,
                     locations = unique(simul_pop$frequencies$location))
  )
  testthat::expect_silent(
    plot_over_time(simul_pop$frequencies, focal_location = 0.5)
  )
  testthat::expect_silent(
    plot_start_end(simul_pop)
  )

  testthat::expect_silent(
    calculate_heterozygosity(simul_pop$population,
                             locations = unique(simul_pop$frequencies$location))
  )

  testthat::expect_silent(
    plot_joyplot_frequencies(simul_pop$frequencies,
                             time_points = c(0, 10, 50))
  )

  testthat::expect_silent(
    calculate_ld(simul_pop$population)
  )
  testthat::expect_silent(
    calculate_marker_frequency(simul_pop, location = 0.5)
  )
})