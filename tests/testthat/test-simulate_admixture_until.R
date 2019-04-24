context("simulate_admixture_until")

test_that("simulate_admixture_until", {

  vx <- simulate_admixture_until(total_runtime = 1000,
                                 pop_size = c(1000, 1000),
                                 morgan = 1,
                                 seed = 42,
                                 generations_between_update = 100,
                                 critical_fst = 0.1,
                                 migration_rate = 0.0)

  fst_2 <- calculate_fst(vx$Population_1,
                         vx$Population_2,
                         sampled_individuals = 100,
                         number_of_markers = 100,
                         random_markers = TRUE)

  testthat::expect_true(vx$FST >= 0.1)
  testthat::expect_true(fst_2 >= 0.1)
  testthat::expect_equal(length(vx$Population_1), 1000)
  testthat::expect_equal(length(vx$Population_2), 1000)
  testthat::expect_true(verify_population(vx$Population_1))
  testthat::expect_true(verify_population(vx$Population_1))
  testthat::expect_true(length(all.equal(vx$Population_1,
                                   vx$Population_2)) > 10)


  vx <- simulate_admixture_until(total_runtime = 1000,
                                 pop_size = c(100, 100),
                                 initial_frequencies = list(c(1,1),
                                                            c(1,1)),
                                 morgan = 1,
                                 seed = 42,
                                 generations_between_update = 100,
                                 critical_fst = 0.1,
                                 migration_rate = 0.01)

  fst_2 <- calculate_fst(vx$Population_1,
                         vx$Population_2,
                         sampled_individuals = 100,
                         number_of_markers = 100,
                         random_markers = TRUE)

  testthat::expect_true(vx$FST >= 0.1)
  testthat::expect_true(fst_2 >= 0.1)
  testthat::expect_equal(length(vx$Population_1), 100)
  testthat::expect_equal(length(vx$Population_2), 100)
  testthat::expect_true(verify_population(vx$Population_1))
  testthat::expect_true(verify_population(vx$Population_1))
  testthat::expect_true(length(all.equal(vx$Population_1,
                                         vx$Population_2)) > 10)




})
