context("simulate_migration")

test_that("simulate_migration", {
  vx <- simulate_admixture_migration(seed = 42,
                                     migration_rate = 0.1)

  testthat::expect_true(verify_population(vx$population_1))
  testthat::expect_true(verify_population(vx$population_2))


  markers <- seq(from = 0.4, to = 0.6, length.out = 100)
  vy <- simulate_admixture_migration(seed = 42,
                                     migration_rate = 0.01,
                                     initial_frequencies = list(c(1,1,0,0),
                                                                c(0,0,1,1)),
                                     total_runtime = 100,
                                     markers = markers)

  testthat::expect_true(verify_population(vy$population_1))
  testthat::expect_true(verify_population(vy$population_2))

  testthat::expect_true( length(markers) == length(unique(vy$frequencies$location)) )


  select_matrix <- matrix(NA, nrow=1, ncol=5)

  s <- 0.1
  select_matrix[1, ] <- c(0.5, 0.5, 0.5+0.5*s, 0.5+s, 0)

  markers <- seq(from = 0.2, to = 0.6, by = 0.001)

  vy <- simulate_admixture_migration(seed = 42,
                                     migration_rate = 0.01,
                                     initial_frequencies = list(c(1,0),
                                                                c(0,1)),
                                     select_matrix = select_matrix,
                                     total_runtime = 100,
                                     markers = markers)

  found <- c()
  for(loc in unique(vy$final_frequency$location)) {
    a <- vy$final_frequency %>%
      filter(location == loc) %>%
      filter(ancestor == 0)
    b <- mean(a$frequency)
    found <- rbind(found, c(loc, b))
  }
  a1 <- subset(found, found[,1] < 0.4)
  a2 <- subset(found, found[,1] == 0.5)

  b1 <- t.test(a1[,2]- 0.5)
  testthat::expect_equal(a2[2], 1)
  testthat::expect_gte(b1$p.value, 0.05)

})
