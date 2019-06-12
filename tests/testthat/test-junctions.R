
test_that("expected_number_junctions", {

  test_expected_junction_number <- function(pop_size,
                                            run_time,
                                            morgan,
                                            replicates)
  {
    cat(pop_size, run_time, morgan, "\n")
    found <- c()
    for (r in 1:replicates) {
      vx <- simulate_admixture(pop_size = pop_size,
                               number_of_founders = 2,
                               total_runtime = run_time,
                               morgan = morgan,
                               seed = r)

      testthat::expect_true(verify_population(vx))

      junct <- calculate_dist_junctions(vx$population)
      found <- c(found, mean(junct))
    }

    require(junctions)
    expected <- junctions::number_of_junctions(N = pop_size,
                                               H_0 = 0.5,
                                               C = morgan,
                                               t = run_time)

    testthat::expect_equal(mean(found), expected, tolerance = 1)
    cat(pop_size, run_time, morgan, mean(found), expected,"\n")
  }

  test_expected_junction_number(pop_size = 100, run_time = 20,
                                morgan = 1, replicates = 20)

  test_expected_junction_number(pop_size = 100, run_time = 100,
                                morgan = 0.5, replicates = 30)

  test_expected_junction_number(pop_size = 100, run_time = 100,
                                morgan = 3, replicates = 30)

  test_expected_junction_number(pop_size = 1000, run_time = 20,
                                morgan = 1, replicates = 30)

  vx <- simulate_admixture(pop_size = 1000, number_of_founders = 2,
                          total_runtime = 5, morgan = 1, seed = 666)

  plot_dist_junctions(vx$population)
})