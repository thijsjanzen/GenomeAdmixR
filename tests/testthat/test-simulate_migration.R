context("simulate_migration")

test_that("simulate_migration", {
  testthat::skip_on_os("solaris")
  message("test simulate_migration")
  vx <- simulate_admixture_migration(migration_rate = 0.1,
                                     total_runtime = 10)

  testthat::expect_true(verify_population(vx$population_1))
  testthat::expect_true(verify_population(vx$population_2))

  testthat::expect_silent(
    vz <- simulate_admixture_migration(input_population_1 = vx$population_1[[1]],
                                     input_population_2 = vx$population_2[[2]],
                                     total_runtime = 10)
  )

  markers <- seq(from = 0.4, to = 0.6, length.out = 100)
  vy <- simulate_admixture_migration(migration_rate = 0.01,
                                     initial_frequencies =
                                       list(c(0.5, 0.5, 0, 0),
                                            c(0, 0, 0.5, 0.5)),
                                     total_runtime = 100,
                                     markers = markers)

  testthat::expect_true(verify_population(vy$population_1))
  testthat::expect_true(verify_population(vy$population_2))

  testthat::expect_true(length(markers) ==
                        length(unique(vy$frequencies$location)))


  if (.Platform$OS.type == "unix") {
      # this seems to not work on windows, I don't know why!

      select_matrix <- matrix(NA, nrow = 1, ncol = 5)

      s <- 0.5
      select_matrix[1, ] <- c(0.5, 0.5, 0.5 + 0.5 * s, 0.5 + s, 0)

      markers <- seq(from = 0.4, to = 0.60, by = 0.01)

      testthat::expect_message(
        vy <- simulate_admixture_migration(migration_rate = 0.01,
                                         initial_frequencies = list(c(1, 0),
                                                                    c(0, 1)),
                                         select_matrix = select_matrix,
                                         total_runtime = 100,
                                         markers = markers)
      )

    found <- c()
    for (loc in unique(vy$final_frequency$location)) {
      a <- subset(vy$final_frequency, vy$final_frequency$location == loc &
                    vy$final_frequency$ancestor == 0)
      b <- mean(a$frequency)
      found <- rbind(found, c(loc, b))
    }
    a2 <- subset(found, found[, 1] == 0.5)

    testthat::expect_equal(a2[2], 1)

    plot_difference_frequencies(vy)
    plot_start_end(vy)

    plot_frequencies(vy$population_1, locations = c(0.3, 0.5, 0.8))
    plot_frequencies(vy$population_2, locations = c(0.3, 0.5, 0.8))
    vv <- plot_joyplot_frequencies(vy$frequencies, time_points = c(0, 10, 50))
  }

  markers <- seq(from = 0.0, to = 1, by = 0.1)

  vy <- simulate_admixture_migration(migration_rate = 0.0,
                                     initial_frequencies =
                                       list(c(0.5, 0.5, 0, 0),
                                            c(0, 0, 0.5, 0.5)),
                                     total_runtime = 100,
                                     markers = markers)

  a1 <- subset(vy$final_frequency, vy$final_frequency$population == 1)
  # this population should only have ancestors 0 and 1

  bv <- c()
  cnt <- 1
  for (i in unique(a1$ancestor)) {
    b1 <- subset(a1, a1$ancestor == i)
    bv[cnt] <- mean(b1$frequency)
    cnt <- cnt + 1
  }

  testthat::expect_equal(bv[3], 0)
  testthat::expect_equal(bv[4], 0)

  a2 <- subset(vy$final_frequency, vy$final_frequency$population == 2)
  bv <- c()
  cnt <- 1
  for (i in unique(a2$ancestor)) {
    b1 <- subset(a2, a2$ancestor == i)
    bv[cnt] <- mean(b1$frequency)
    cnt <- cnt + 1
  }
  testthat::expect_equal(bv[1], 0)
  testthat::expect_equal(bv[2], 0)

  vy <- simulate_admixture_migration(migration_rate = 0.0,
                                     initial_frequencies =
                                       list(c(0.5, 0.5, 0, 0),
                                            c(0, 0, 0.5, 0.5)),
                                     total_runtime = 100,
                                     markers = 0.5,
                                     track_junctions = TRUE)
})

test_that("simulate_migration no seed", {
  testthat::skip_on_os("solaris")
  message("test simulate_migration no seed")
   # no markers:
  testthat::expect_warning(
    vy <- simulate_admixture_migration(migration_rate = 0.0,
                                       initial_frequencies =
                                         list(c(1, 1, 0, 0),
                                              c(0, 0, 1, 1)),
                                       total_runtime = 100,
                                       track_junctions = TRUE)
  )
})
