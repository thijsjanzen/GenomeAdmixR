context("tajima")

testthat::test_that("tajima", {
  pop <- simulate_admixture(pop_size = 100,
                            number_of_founders = 2,
                            seed = 666,
                            total_runtime = 10)

  t <- 10

  found <- isoSIM::calculate_tajima_d(pop$population)$D

  while(t < 1000) {
    pop <- simulate_admixture(pop$population,
                              pop_size = 100,
                              seed = t,
                              total_runtime = 10,
                              progress_bar = FALSE)
    found <- c(found,
               isoSIM::calculate_tajima_d(pop$population)$D
    )
    t <- t + 10
   # cat(t,"\n")
  }

  testthat::expect_true(mean(found, na.rm = T) < 2)
})

testthat::test_that("tajima abuse", {
  pop <- simulate_admixture(pop_size = 100,
                            number_of_founders = 2,
                            seed = 666,
                            total_runtime = 10)

  testthat::expect_error(  calculate_tajima_d(pop$initial_frequency) )
})


testthat::test_that("tajima", {

  if(1 == 2) {

  pop_size <- 1000
  found <- c()

  for(r in 1:10) {

  pop <- c()
  if(r ) {
    cat("no selection\n")
      pop <- simulate_admixture(pop_size = pop_size,
                            number_of_founders = 4,
                            seed = 666+r,
                            total_runtime = 1000)
  } else {
    cat("selection\n")
      select_matrix <- matrix(nrow=1,ncol=5)
      select_matrix[1,] <- c(0.5, 1, 1.1, 1.2, 0)
      pop <- simulate_admixture(pop_size = pop_size,
                                number_of_founders = 4,
                                seed = 666+r,
                                select_matrix = select_matrix,
                                total_runtime = 1000)
  }
  pop <- pop$population

  markers <- seq(0,1,length.out = 100)

  number_of_sampled_individuals <- pop_size

  D1 <- calculate_tajima_d(pop,
                                   markers,
                                   number_of_sampled_individuals)

  # now prepare data for comparison

  indices_sampled_individuals <- sample(1:pop_size,
                                        number_of_sampled_individuals,
                                        replace = F)

  loci_matrix <- matrix(ncol = length(markers),
                        nrow = 2*number_of_sampled_individuals)

  for(m in seq_along(markers)) {
    for(indiv in seq_along(indices_sampled_individuals)) {

      indiv_start <- (indiv*2)-1

      loci_matrix[indiv_start, m] <- findtype(
        pop[[indices_sampled_individuals[indiv]]]$chromosome1,
        markers[m])
      loci_matrix[indiv_start + 1, m] <- findtype(
        pop[[indices_sampled_individuals[indiv]]]$chromosome2,
        markers[m])
    }
  }

  sim_data <- list()
  for(i in seq_along(loci_matrix[,1])) {
    x <- loci_matrix[i,]
    x <- x[!is.na(x)]
    x[x==0] <- "t"
    x[x==1] <- "c"
    x[x==2] <- "g"
    x[x==3] <- "a"
    sim_data[[i]] <- x
  }

  require(strataG)
  D2 <- strataG::tajimasD(sim_data)

 # testthat::expect_equal(D1$D, D2[1], tolerance = 0.1)
  cat(r, D1$D, D2, "\n")
  found <- rbind(found, c(D1$D, D2[1]))
  }
  }


})



