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
                            morgan = morgan)


  testthat::expect_true(verify_population(pop))

testthat::expect_silent(
  females <- create_iso_female(pop, n = 5, run_time = 3000)
)
testthat::expect_silent(
  females <- create_iso_female(pop, n = 1, run_time = 3000)
)

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
                             morgan = morgan)

  pop2 <- simulate_admixture(pop_size = pop_size,
                             number_of_founders = number_of_founders,
                             total_runtime = run_time,
                             morgan = morgan)
  pop2 <- increase_ancestor(pop2, number_of_founders)

  testthat::expect_true(verify_population(pop1))
  testthat::expect_true(verify_population(pop2))


  testthat::expect_silent(
      female_1 <- create_iso_female(pop1, n = 1,
                                run_time = 20000)
  )
  testthat::expect_silent(
   female_2 <- create_iso_female(pop2, n = 1,
                                run_time = 20000)
)
  testthat::expect_true(verify_individual(female_1[[1]]))
  testthat::expect_true(verify_individual(female_2[[1]]))

testthat::expect_silent(
  females <- create_iso_female(pop1, n = 2, run_time = 2000)
)
  pop_size = 100
  vy <- simulate_admixture(input_population = females,
                           pop_size = pop_size,
                           total_runtime = 200,
                           morgan = morgan)

  testthat::expect_equal(length(vy$population), pop_size)
  testthat::expect_true(verify_population(vy))

  vy <- simulate_admixture(list(female_1[[1]],  female_2[[1]]),
                           pop_size,
                           2000,
                           morgan)

  testthat::expect_equal(length(vy$population), pop_size)
  testthat::expect_true(verify_population(vy))

  testthat::expect_silent(
    plot_chromosome(female_1[[1]]$chromosome1, 0, 1)
  )
})

test_that("cpp classes", {
  a <- matrix(c(0.1, 1, 2, 2), nrow = 2)
  b <- matrix(c(0, 1, 1, -1), nrow = 2)
  indiv <- list(chromosome1 = a, chromosome2 = a)
  class(indiv) <- "individual"

  # chromosome 1
  testthat::expect_warning(v <- verify_individual(indiv),
                           "Chromosome doesn't start at 0")

  indiv <- list(chromosome1 = b, chromosome2 = a)
  class(indiv) <- "individual"

  # chromosome 2
  testthat::expect_warning(v <- verify_individual(indiv),
                          "Chromosome doesn't start at 0")

  a <- matrix(c(0.0, 1, 2, 2), nrow = 2)
  b <- matrix(c(0, 1, 1, -1), nrow = 2)
  indiv <- list(chromosome1 = b, chromosome2 = a)
  class(indiv) <- "individual"
  testthat::expect_warning(v <- verify_individual(indiv),
                          "Chromosome doesn't end with -1")
  indiv$chromosome2 <-  indiv$chromosome1
  indiv$chromosome1 <- a
  testthat::expect_warning(v <- verify_individual(indiv),
                          "Chromosome doesn't end with -1")


  a <- matrix(c(0.0, 1, 0.5, 29192875037,  1, -1), ncol = 2)

  indiv$chromosome1 <- a
  testthat::expect_warning(v <- verify_individual(indiv),
                          "Memory error recorded in chromosome")

  a <- matrix(c(0.0, 1, 0.5, -92875037,  1, -1), ncol = 2)
  indiv$chromosome2 <- a
  indiv$chromosome1 <- b
  testthat::expect_warning(v <- verify_individual(indiv),
                          "Memory error recorded in chromosome")
})
