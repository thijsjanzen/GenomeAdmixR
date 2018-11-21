calculate_tajima_d <- function(pop,
                               markers = seq(1e-6,1-1e-6,length.out = 100),
                               number_of_sampled_individuals = 10) {

  pop <- check_input_pop(pop)

  pop_size <- length(pop)
  indices_sampled_individuals <- sample(1:pop_size,
                                        number_of_sampled_individuals,
                                        replace = FALSE)

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

  mutation_rate <- 0
  expected_theta <- 4*pop_size * mutation_rate

  n <- 2*number_of_sampled_individuals
  tmp <- 1:(n - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (n + 1)/(3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)

  pi <- 0 #average number of pairwise differences
  cnt <- 0
  for(i in seq_along(loci_matrix[,1])) {
    for(j in 1:i) {
      if(i != j) {
        loci_1 <- loci_matrix[i,]
        loci_2 <- loci_matrix[j,]
        difference <- loci_1 - loci_2
        num_diff <- length(difference) - length(which(difference == 0))
        pi <- pi + num_diff
        cnt <- cnt + 1
        #    cat(i,j,"\n")
      }
    }
  }

  if(choose(2*number_of_sampled_individuals, 2) != cnt) {
    stop("Too many pairwise comparisons!\n")
  }

  pi <- pi * 1.0 / (choose(2*number_of_sampled_individuals, 2))

  S <- 0 # number of segregating sites
  for(i in seq_along(loci_matrix[1,])) {
    a <- as.numeric(loci_matrix[,i])
    if(length(unique(a)) > 1) {
      S <- S + 1
    }
  }

  D <- (pi - S / a1) / sqrt(e1*S + e2 * S * (S - 1))

  Dmin <- (2/n - 1/a1)/sqrt(e2)
  Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
  tmp1 <- 1 + Dmin * Dmax
  tmp2 <- Dmax - Dmin
  a <- -tmp1 * Dmax/tmp2
  b <- tmp1 * Dmin/tmp2
  p <- pbeta((D - Dmin)/tmp2, b, a)
  if(!is.nan(p)) {
    p <- if (p < 0.5) {
      2 * p } else {2 * (1 - p)}
  }

  theta_hat <- S / a1

  return(list("D" = D,
              "Pi" = pi,
              "S" = S,
              "theta_hat[Estimated_from_S]" = theta_hat,
              "p value" = p))
}