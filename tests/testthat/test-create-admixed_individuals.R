context("create_admixed_individuals")

testthat::test_that("create_admixed_individuals", {
  admixed_pop <- create_admixed_individuals(num_individuals = 10,
                                            population_size = 100,
                                            number_of_founders = 2,
                                            size_in_morgan = 1)

  testthat::expect_true( is(admixed_pop$population, "population") )

  expected_num_junctions <- junctions::number_of_junctions(N = 100,
                                                           R = Inf,
                                                           H_0 = 0.5,
                                                           C = 1,
                                                           t = Inf)
  found <- c()
  for(i in seq_along(admixed_pop$population)) {
    found <- c(found, length(admixed_pop$population[[i]]$chromosome1[,2]) - 2)
    found <- c(found, length(admixed_pop$population[[i]]$chromosome2[,2]) - 2)
  }
  avg_junctions <- mean(found)
  testthat::expect_equal(expected_num_junctions, avg_junctions, tolerance = 5)
})
