context("general usage")


test_that("general usage", {

  data("dgrp2.3R.5k.data")

  mks <- sample(dgrp2.3R.5k.data$markers, size = 100,
               replace = FALSE, prob = NULL)


  testthat::expect_silent(
    simulated_pop <- simulate_admixture_data(input_data = dgrp2.3R.5k.data,
                                           pop_size = 1000,
                                           total_runtime = 10,
                                           morgan = 1,
                                           num_threads = -1,
                                           markers = mks)
  )

  ### The list option is working now with genomeadmixr_data type:
  testthat::expect_message(
  simulated_pop_2 <- simulate_admixture_data(input_data =
                                               list(dgrp2.3R.5k.data,
                                                    dgrp2.3R.5k.data),
                                             pop_size = 1000,
                                             total_runtime = 10,
                                             morgan = 1,
                                             markers = mks),
     "found multiple input populations"
  )

  testthat::expect_message(
  simulated_pop_2 <- simulate_admixture_data(input_data =
                                               list(simulated_pop,
                                                    simulated_pop),
                                             pop_size = 1000,
                                             total_runtime = 10,
                                             verbose = TRUE,
                                             num_threads = -1,
                                             morgan = 1,
                                             markers = mks)
  )

  calculate_marker_frequency(simulated_pop, location = mks[50])
  subset(simulated_pop$frequencies, time == 10 & location == mks[50])

  selection_matrix <- matrix(nrow = 1, ncol = 5)
  selection_matrix[1, ] <- c(mks[50], 0.4, 0.7, 1.0, "t")


  testthat::expect_silent(
      iso_100 <- create_iso_female_data(input_data = dgrp2.3R.5k.data,
                                   inbreeding_pop_size = 100,
                                   n = 20,
                                   morgan = 1,
                                   run_time = 20)
  )

  testthat::expect_message(
    iso_100 <- create_iso_female_data(input_data = simulated_pop,
                                   inbreeding_pop_size = 100,
                                   n = 20,
                                   morgan = 1,
                                   num_threads = 2,
                                   run_time = 20)
  )

  testthat::expect_message(
    selected_pop <- simulate_admixture_data(input_data = simulated_pop,
                                          pop_size = 1000,
                                          total_runtime = 16,
                                          morgan = 1,
                                          markers = mks,
                                          select_matrix = selection_matrix)
  )

  testthat::expect_message(
    two_pops <- simulate_admixture_migration_data(input_data_population_1 =
                                                 simulated_pop,
                                               input_data_population_2 =
                                                 simulated_pop,
                                               pop_size = 100,
                                               total_runtime = 10,
                                               morgan = 1,
                                               migration_rate = 0,
                                               markers = mks,
                                               stop_at_critical_fst = TRUE,
                                               critical_fst = 0.05,
                                               generations_between_update = 100,
                                               num_threads = 4)
  )
})

test_that("isofemale usage", {
  data("dgrp2.3R.5k.data")

  mks = sample(dgrp2.3R.5k.data$markers, size = 300, replace = FALSE, prob = NULL)

  testthat::expect_silent(
  simulated_pop <- simulate_admixture_data(input_data = dgrp2.3R.5k.data,
                                           pop_size = 1000,
                                           total_runtime = 10,
                                           morgan = 1,
                                           markers = mks)
  )

  testthat::expect_message(
  isos <- create_iso_female_data(input_data = simulated_pop,
                                 n = 20,
                                 inbreeding_pop_size = 100,
                                 run_time = 50)
  )

  testthat::expect_message(
    simulated_pop <- simulate_admixture_data(input_data = isos[1],
                                           pop_size = 1000,
                                           total_runtime = 10,
                                           morgan = 1,
                                           markers = mks)
  )
})