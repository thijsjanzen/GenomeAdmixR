# Written by Claude (Anthropic): regression tests for the bugs fixed in 2.1.14.

make_f1_data <- function(n_markers, n_ind = 50) {
  # every individual carries one chromosome of all 'a' and one of all 'c'
  g <- matrix(0, 2 * n_ind, n_markers)
  g[seq(1, 2 * n_ind, 2), ] <- 1
  g[seq(2, 2 * n_ind, 2), ] <- 2
  d <- list(genomes = g, markers = seq_len(n_markers) * 1000)
  class(d) <- "genomeadmixr_data"
  return(d)
}

test_that("tracked junctions are counted correctly", {
  vx <- simulate_admixture(module = ancestry_module(morgan = 1,
                                                   track_junctions = TRUE),
                           pop_size = 50, total_runtime = 5)
  testthat::expect_true(all(vx$junctions >= 0))
  testthat::expect_true(all(vx$junctions < 10))

  vy <- simulate_admixture(module = ancestry_module(morgan = 1e-12,
                                                   track_junctions = TRUE),
                           pop_size = 50, total_runtime = 3)
  testthat::expect_equal(vy$junctions, rep(0, length(vy$junctions)))
})

test_that("selection acts on chromosomes without junctions", {
  # without recombination chromosomes remain a single ancestry block
  sel <- matrix(c(0.5, 1, 2, 4, 0), nrow = 1)
  f <- replicate(3, {
    s <- simulate_admixture(module = ancestry_module(number_of_founders = 2,
                                                     morgan = 1e-6,
                                                     markers = 0.5),
                            select_matrix = sel, pop_size = 200,
                            total_runtime = 30)
    ff <- s$final_frequency
    ff$frequency[ff$ancestor == 0]
  })
  testthat::expect_gt(mean(f), 0.9)
})

test_that("selection site is in Morgan", {
  sel <- matrix(c(1.5, 1, 2, 4, 0), nrow = 1)
  s <- simulate_admixture(module = ancestry_module(number_of_founders = 2,
                                                   morgan = 2,
                                                   markers = c(1.5)),
                          select_matrix = sel, pop_size = 200,
                          total_runtime = 30)
  ff <- s$final_frequency
  testthat::expect_gt(ff$frequency[ff$ancestor == 0], 0.9)
})

test_that("vcf genotypes map to the right haplotypes", {
  # three individuals: hom ref, het, hom alt
  m <- vcf_to_matrix_cpp(matrix(c(1, 2, 3), ncol = 1), 1, 2)
  testthat::expect_equal(as.vector(m), c(1, 1, 1, 2, 2, 2))
})

test_that("default recombination map is uniform along the chromosome", {
  d <- make_f1_data(200)
  s <- simulate_admixture(module = sequence_module(molecular_data = d,
                                                   morgan = 1),
                          pop_size = 2000, total_runtime = 1)
  rel_pos <- unlist(lapply(s$population, function(ind) {
    lapply(list(ind$chromosome1, ind$chromosome2), function(ch) {
      which(diff(ch[, 2]) != 0) / nrow(ch)
    })
  }))
  testthat::expect_equal(mean(rel_pos), 0.5, tolerance = 0.05)
})

test_that("recombination map is not reused across simulations", {
  for (k in 1:2) {
    simulate_admixture(module = sequence_module(molecular_data =
                                                  make_f1_data(2000),
                                                morgan = 3),
                       pop_size = 500, total_runtime = 2, num_threads = 2)
    s <- simulate_admixture(module = sequence_module(molecular_data =
                                                       make_f1_data(10),
                                                     morgan = 3),
                            pop_size = 500, total_runtime = 2,
                            num_threads = 2)
  }
  testthat::expect_true(all(vapply(s$population, function(i) {
    nrow(i$chromosome1) == 10 && nrow(i$chromosome2) == 10
  }, TRUE)))
})

test_that("allele frequencies do not depend on the other locations", {
  pop <- simulate_admixture(module = ancestry_module(number_of_founders = 2),
                            pop_size = 100, total_runtime = 1)$population
  a <- calculate_allele_frequencies(pop, locations = c(0.1, 0.5))
  b <- calculate_allele_frequencies(pop, locations = c(0.1, 0.5, 1))
  testthat::expect_equal(a$frequency, b$frequency[b$location < 1])
})

test_that("multiple selected sites in sequence data", {
  d <- make_f1_data(20)
  sel <- matrix(c(5000, 1, 1, 1, "a",
                  9000, 1, 1, 1, "a"), nrow = 2, byrow = TRUE)
  testthat::expect_silent(
    suppressMessages(
      simulate_admixture(module = sequence_module(molecular_data = d,
                                                  morgan = 1),
                         select_matrix = sel, pop_size = 50,
                         total_runtime = 2)
    )
  )
})

test_that("early fixation does not leave empty frequency rows", {
  s <- simulate_admixture(module = ancestry_module(number_of_founders = 2,
                                                   morgan = 1e-6,
                                                   markers = 0.5),
                          pop_size = 6, total_runtime = 500)
  testthat::expect_false(any(s$frequencies$time == 0 &
                               s$frequencies$location == 0))
  testthat::expect_equal(s$final_frequency$time[[1]],
                         max(s$frequencies$time) + 1)
})

test_that("simulations are reproducible with set.seed", {
  run <- function(threads) {
    set.seed(42)
    simulate_admixture(module = ancestry_module(), pop_size = 100,
                       total_runtime = 10, num_threads = threads)$population
  }
  testthat::expect_identical(run(2), run(2))
  testthat::expect_identical(run(1), run(2))
})

test_that("input_population is used", {
  p <- simulate_admixture(module = ancestry_module(number_of_founders = 2,
                                                   initial_frequencies =
                                                     c(1, 0)),
                          pop_size = 50, total_runtime = 1)$population
  s <- simulate_admixture(module = ancestry_module(input_population = p,
                                                   markers = 0.5),
                          pop_size = 50, total_runtime = 1)
  f <- s$initial_frequency
  testthat::expect_equal(f$frequency[f$ancestor == 0], 1)
})

test_that("calculate_ld with few individuals", {
  pop <- simulate_admixture(pop_size = 5, total_runtime = 2)$population
  testthat::expect_silent(calculate_ld(pop, sampled_individuals = 10,
                                       markers = 5))
})
