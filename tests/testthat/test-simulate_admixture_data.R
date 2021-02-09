context("test simulate admixture data")

test_that("simulate_admixture_data", {
  testthat::skip_on_os("solaris")
  message("test simulate_admixture_data")

  num_markers <- 100
  num_indiv <- 100
  chosen_markers <- 1:num_markers

  fake_input_data1 <- create_artificial_genomeadmixr_data(
    number_of_individuals = num_indiv,
    marker_locations = chosen_markers,
    used_nucleotides = 1:2
  )

  fake_input_data2 <- create_artificial_genomeadmixr_data(
    number_of_individuals = num_indiv,
    marker_locations = chosen_markers,
    used_nucleotides = 3:4
  )

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

test_that("simulate_admixture_data_mutation", {
  testthat::skip_on_os("solaris")
  message("test simulate_admixture_data_mutation")

  num_markers <- 100
  num_indiv <- 100
  chosen_markers <- 1:num_markers

  fake_input_data1 <- create_artificial_genomeadmixr_data(
    number_of_individuals = num_indiv,
    marker_locations = chosen_markers,
    used_nucleotides = 1:2
  )

  fake_input_data2 <- create_artificial_genomeadmixr_data(
    number_of_individuals = num_indiv,
    marker_locations = chosen_markers,
    used_nucleotides = 3:4
  )

  combined_data <- combine_input_data(input_data_list = list(fake_input_data1,
                                                             fake_input_data2),
                                      frequencies = c(0.5, 0.5),
                                      pop_size = 1000)

  sub_matrix <- matrix(0.25, nrow = 4, ncol = 4)

  testthat::expect_message(
    simul_pop <- simulate_admixture_data(input_data = combined_data,
                                       pop_size = 100,
                                       total_runtime = 100,
                                       markers = chosen_markers,
                                       morgan = 1,
                                       verbose = FALSE,
                                       mutation_rate = 0.1,
                                       substitution_matrix = sub_matrix)
  )

  a1 <- simul_pop$initial_frequency
  a2 <- simul_pop$final_frequency

  aa <- a1 %>%
          dplyr::group_by(ancestor) %>%
    dplyr::summarise("mean_freq" = mean(frequency))

  bb <- a2 %>%
    dplyr::group_by(ancestor) %>%
    dplyr::summarise("mean_freq" = mean(frequency))

  testthat::expect_equal(mean(bb$mean_freq[2:5]), 0.25)

  bases <- c("a", "c", "t", "g")
  for (i in 1:4) {
    sub_matrix <- matrix(0, nrow = 4, ncol = 4)
    sub_matrix[, i] <- 1

    testthat::expect_message(
      testthat::expect_warning(
      simul_pop <- simulate_admixture_data(input_data = combined_data,
                                           pop_size = 100,
                                           total_runtime = 100,
                                           markers = chosen_markers,
                                           morgan = 1,
                                           verbose = FALSE,
                                           mutation_rate = 1.0,
                                           substitution_matrix = sub_matrix)
      )
    )

    bb <- simul_pop$final_frequency %>%
      dplyr::group_by(ancestor) %>%
      dplyr::summarise("mean_freq" = mean(frequency))

    highest_base <- bb$ancestor[which.max(bb$mean_freq)]
    testthat::expect_equal(bases[i], highest_base)
  }
})
