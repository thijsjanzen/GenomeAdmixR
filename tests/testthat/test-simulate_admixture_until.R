context("simulate_admixture_until")

test_that("simulate_admixture_until", {
  testthat::skip_on_os("solaris")
 message("test simulate_admixture_until")
    vx <- simulate_admixture_migration(total_runtime = 1000,
                                     pop_size = c(100, 100),
                                     initial_frequencies = list(c(0.5, 0.5),
                                                                c(0.5, 0.5)),
                                     morgan = 1,
                                     stop_at_critical_fst = TRUE,
                                     generations_between_update = 10,
                                     critical_fst = 0.2,
                                     migration_rate = 0.001)

  fst_2 <- calculate_fst(vx$population_1,
                         vx$population_2,
                         sampled_individuals = 100,
                         number_of_markers = 100,
                         random_markers = TRUE)

  testthat::expect_true(vx$FST >= 0.05)
  testthat::expect_true(fst_2 >= 0.05)


  testthat::expect_equal(length(vx$population_1), 100)
  testthat::expect_equal(length(vx$population_2), 100)
  testthat::expect_true(verify_population(vx$population_1))
  testthat::expect_true(verify_population(vx$population_2))
  testthat::expect_true(length(all.equal(vx$population_1,
                                         vx$population_2)) > 10)
})

test_that("simulate_admixture_until_data", {
  testthat::skip_on_os("solaris")
  message("test simulate_admixture_until_data")

  num_markers <- 100
  num_indiv <- 100
  chosen_markers <- 1:num_markers

  fake_input_data1 <- list()
  fake_input_data1$genomes <- matrix(data = 1,
                                     nrow = num_indiv,
                                     ncol = num_markers)


  fake_input_data1$markers <- chosen_markers

  fake_input_data2 <- list()
  fake_input_data2$genomes <- matrix(data = 2,
                                     nrow = num_indiv,
                                     ncol = num_markers)
  fake_input_data2$markers <- chosen_markers

  class(fake_input_data1) <- "genomeadmixr_data"
  class(fake_input_data2) <- "genomeadmixr_data"

  vx <- simulate_admixture_migration_data(
    input_data_population_1 = fake_input_data1,
    input_data_population_2 = fake_input_data2,
    pop_size = c(100, 100),
    total_runtime = 100,
    markers = chosen_markers,
    morgan = 1,
    migration_rate = 0.001,
    critical_fst = 0.2,
    generations_between_update = 10,
    verbose = FALSE)

  fst_2 <- calculate_fst(vx$population_1,
                         vx$population_2,
                         sampled_individuals = 10,
                         number_of_markers = 100,
                         random_markers = TRUE)

  testthat::expect_true(vx$FST >= 0.05)
  testthat::expect_true(fst_2 >= 0.05)


  testthat::expect_equal(length(vx$population_1), 100)
  testthat::expect_equal(length(vx$population_2), 100)
  testthat::expect_true(length(all.equal(vx$population_1,
                                         vx$population_2)) > 10)
})

