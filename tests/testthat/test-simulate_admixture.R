context("simulate_Admixture")


test_that("simulate admixture use, markers", {
  pop <- simulate_admixture(pop_size = 1000,
                            number_of_founders = 2,
                            total_runtime = 3,
                            markers = seq(0, 1, length.out = 1000),
                            morgan = 1)

  a <- subset(pop$final_frequency, pop$final_frequency$location < 0.5)
  b <-  subset(pop$final_frequency, pop$final_frequency$location > 0.5)
  testthat::expect_equal(dim(a), dim(b))

  testthat::expect_silent(
    GenomeAdmixR::plot_difference_frequencies(pop)
  )


  pop <- simulate_admixture(pop_size = 1000,
                            number_of_founders = 2,
                            total_runtime = 3,
                            markers = seq(0, 3, length.out = 1000),
                            morgan = 3)

  a <- subset(pop$final_frequency, pop$final_frequency$location < 1.5)
  b <-  subset(pop$final_frequency, pop$final_frequency$location > 1.5)
  testthat::expect_equal(dim(a), dim(b))

  GenomeAdmixR::plot_difference_frequencies(pop)

  pop <- simulate_admixture(pop_size = 1000,
                            number_of_founders = 2,
                            total_runtime = 3,
                            markers = seq(0, 5, length.out = 1000),
                            morgan = 5)

  a <- subset(pop$final_frequency, pop$final_frequency$location < 2.5)
  b <-  subset(pop$final_frequency, pop$final_frequency$location > 2.5)
  testthat::expect_equal(dim(a), dim(b))
  GenomeAdmixR::plot_difference_frequencies(pop)
})