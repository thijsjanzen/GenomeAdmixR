context("isoFemale creation")

test_that("create_isofemale", {

  pop_size <- 100
  number_of_founders <- 2
  run_time <- 100
  morgan <- 1
  write_to_file <- FALSE

  pop <- simulate_admixture(pop_size = pop_size,
                            number_of_founders = number_of_founders,
                            total_runtime = run_time,
                            morgan = morgan,
                            seed = 42)


  testthat::expect_true(verify_population(pop))

  females <- create_iso_female(pop, n = 5, run_time = 3000)

  females <- create_iso_female(pop, n = 1, run_time = 3000)


  testthat::expect_equal(length(females), 1)
})

test_that("create_population_from_isofemales", {

  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.5

  pop1 <- simulate_admixture(pop_size = pop_size,
                             number_of_founders = number_of_founders,
                             total_runtime = run_time,
                             morgan = morgan,
                             seed = 42)

  pop2 <- simulate_admixture(pop_size = pop_size,
                             number_of_founders = number_of_founders,
                             total_runtime = run_time,
                             morgan = morgan,
                             seed = 24)
  pop2 <- increase_ancestor(pop2, number_of_founders)

  testthat::expect_true(verify_population(pop1))
  testthat::expect_true(verify_population(pop2))

  female_1 <- create_iso_female(pop1, n = 1,
                                run_time = 2000, seed = 1)
  female_2 <- create_iso_female(pop2, n = 1,
                                run_time = 2000, seed = 2)



  testthat::expect_true(verify_individual(female_1[[1]]))
  testthat::expect_true(verify_individual(female_2[[1]]))


  females <- create_iso_female(pop1, n = 2, run_time = 2000)

  vy <- simulate_admixture(females,
                           pop_size, 2000,
                           morgan,
                           seed = 666)

  testthat::expect_equal(length(vy$population), pop_size)
  testthat::expect_true(verify_population(vy))

  vy <- simulate_admixture(list(female_1[[1]],  female_2[[1]]),
                           pop_size,
                           2000,
                           morgan,
                           seed = 666)

  testthat::expect_equal(length(vy$population), pop_size)
  testthat::expect_true(verify_population(vy))

  plot_chromosome(female_1[[1]]$chromosome1, 0, 1)
})

test_that("cpp classes", {
  test_fish_functions()

  a <- matrix(c(0.1, 1, 2, 2), nrow = 2)
  b <- matrix(c(0, 1, 1, -1), nrow = 2)
  indiv <- list(chromosome1 = a, chromosome2 = a)
  class(indiv) <- "individual"

  # chromosome 1
  testthat::expect_output(v <- verify_individual(indiv),
                           "Chromosome doesn't start at 0")

  indiv <- list(chromosome1 = b, chromosome2 = a)
  class(indiv) <- "individual"

  # chromosome 2
  testthat::expect_output(v <- verify_individual(indiv),
                          "Chromosome doesn't start at 0")

  a <- matrix(c(0.0, 1, 2, 2), nrow = 2)
  b <- matrix(c(0, 1, 1, -1), nrow = 2)
  indiv <- list(chromosome1 = b, chromosome2 = a)
  class(indiv) <- "individual"
  testthat::expect_output(v <- verify_individual(indiv),
                          "Chromosome doesn't end with -1")
  indiv$chromosome2 <-  indiv$chromosome1
  indiv$chromosome1 <- a
  testthat::expect_output(v <- verify_individual(indiv),
                          "Chromosome doesn't end with -1")


  a <- matrix(c(0.0, 1, 0.5, 29192875037,  1, -1), ncol = 2)

  indiv$chromosome1 <- a
  testthat::expect_output(v <- verify_individual(indiv),
                          "Memory error recorded in chromosome")

  a <- matrix(c(0.0, 1, 0.5, -92875037,  1, -1), ncol = 2)
  indiv$chromosome2 <- a
  indiv$chromosome1 <- b
  testthat::expect_output(v <- verify_individual(indiv),
                          "Memory error recorded in chromosome")
})
