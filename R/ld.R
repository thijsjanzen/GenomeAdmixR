#nolint start
calculate_average_LD <- function(alleles_pos_1, alleles_pos_2) {

  all_alleles <- c(as.vector(alleles_pos_1), as.vector(alleles_pos_2))
  all_alleles <- sort(unique(all_alleles))

  LD <- 0
  r_squared <- 0
  for (i in all_alleles) {
    for (j in all_alleles) {
      p_A_i <- length(which(alleles_pos_1 == i)) / length(alleles_pos_1)
      p_B_j <- length(which(alleles_pos_2 == j)) / length(alleles_pos_2)

      countAB <- 0
      for (a in seq_along(alleles_pos_1[, 1])) {
        if ( (alleles_pos_1[a, 1]) == i && (alleles_pos_2[a, 1] == j)) {
          countAB <- countAB + 1
        }

        if ( (alleles_pos_1[a, 2]) == i && (alleles_pos_2[a, 2] == j)) {
          countAB <- countAB + 1
        }
      }

      p_A_i_B_j <- countAB / (length(alleles_pos_1[, 1]) +
                              length(alleles_pos_1[, 2]))

      if (is.nan(p_A_i_B_j)) p_A_i_B_j <- 0

      D_i_j <- p_A_i_B_j - p_A_i * p_B_j

      D_i_j_max <- 1

      if (D_i_j < 0) {
        D_i_j_max <- min(p_A_i * p_B_j, (1 - p_A_i) * (1 - p_B_j))
      } else {
        D_i_j_max <- min(p_A_i * (1 - p_B_j), (1 - p_A_i) * p_B_j)
      }

      if (D_i_j_max > 0) {
        LD <- LD + p_A_i * p_B_j * abs(D_i_j / D_i_j_max)
        r_i_j <- (D_i_j ^ 2) / (p_A_i * (1 - p_A_i) * p_B_j * (1 - p_B_j))
        r_squared <- r_squared + p_A_i * p_B_j * r_i_j
      }
    }
  }

  return(list("LD" = LD,
              "r_sq" = r_squared)
         )
}

calculate_LD <- function(pop,
                         sampled_individuals = 10,
                         number_of_markers = 100,
                         random_markers = TRUE) {

  pop <- check_input_pop(pop)

  all_loci <- matrix(nrow = length(pop), ncol = 2 * number_of_markers, 0)

  markers <- seq(1e-9, 1 - (1e-9), length.out = number_of_markers)
  if (random_markers) {
    markers <- create_random_markers(number_of_markers)
  }

  for (x in seq_along(markers)) {
    focal_marker <- markers[x]
    for (i in seq_along(pop)) {
      allele_1 <- 1 + findtype(pop[[i]]$chromosome1, focal_marker)
      allele_2 <- 1 + findtype(pop[[i]]$chromosome2, focal_marker)

      index <- (x - 1) * 2 + 1

      all_loci[i, index]     <- as.numeric(allele_1)
      all_loci[i, index + 1] <- as.numeric(allele_2)
    }
  }

  LD_matrix   <- matrix(nrow = length(markers), ncol = length(markers), NA)
  rsq_matrix  <- matrix(nrow = length(markers), ncol = length(markers), NA)
  dist_matrix <- matrix(nrow = length(markers), ncol = length(markers), NA)

  for (x in seq_along(markers)) {
    for (y in seq_len(x)) {
      if (x != y) {
        index1 <- c( (x - 1) * 2 + 1, (x - 1) * 2 + 2)
        index2 <- c( (y - 1) * 2 + 1, (y - 1) * 2 + 2)
        g1 <- all_loci[, index1]
        g2 <- all_loci[, index2]

        ld <- calculate_average_LD(g1, g2)
        LD_matrix[x, y] <- ld$LD
        rsq_matrix[x, y] <- ld$r_sq
        gen_dist <- abs(markers[x] - markers[y])
        dist_matrix[x, y] <- gen_dist
      }
    }
  }

  return(list("LD_matrix" = LD_matrix,
              "rsq_matrix" = rsq_matrix,
              "dist_matrix" = dist_matrix))
}
#nolint end