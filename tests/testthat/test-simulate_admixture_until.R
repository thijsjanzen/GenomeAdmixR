context("simulate_admixture_until")

test_that("simulate_admixture_until", {
 #skip("untilll")
  testthat::expect_message(
    vx <- simulate_admixture_migration(total_runtime = 1000,
                                     pop_size = c(100, 100),
                                     initial_frequencies = list(c(0.5, 0.5),
                                                                c(0.5, 0.5)),
                                     morgan = 1,
                                     stop_at_critical_fst = TRUE,
                                     generations_between_update = 10,
                                     critical_fst = 0.2,
                                     migration_rate = 0.001)
  )

  fst_2 <- calculate_fst(vx$Population_1,
                         vx$Population_2,
                         sampled_individuals = 100,
                         number_of_markers = 100,
                         random_markers = TRUE)

  testthat::expect_true(vx$FST >= 0.05)
  testthat::expect_true(fst_2 >= 0.05)


  testthat::expect_equal(length(vx$Population_1), 100)
  testthat::expect_equal(length(vx$Population_2), 100)
  testthat::expect_true(verify_population(vx$Population_1))
  testthat::expect_true(verify_population(vx$Population_1))
  testthat::expect_true(length(all.equal(vx$Population_1,
                                         vx$Population_2)) > 10)
})
