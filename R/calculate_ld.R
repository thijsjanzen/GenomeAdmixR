#' Calculate linkage disequilibrium statistics
#' This function calculates two matrices, once containing all pairwise
#' linkage disequilibrium (ld) values, and one matrix containing all pairwise r
#' statistics
#' @param pop focal population
#' @param sampled_individuals Number of individuals randomly sampled to
#' calculate the LD matrices
#' @param number_of_markers Number of markers used to calculate the ld matrices
#' @param random_markers If TRUE, markers are randomly spaced along the
#' chromosome, if FALSE, markers are equidistantly spaced along the chromosome.
#' @return An object containing two items:
#' \item{ld_matrix}{
#'   Pairwise ld statistics for all markers
#' }
#' \item{rsq_matrix}{
#'   Pairwise rsq statistics for all markers
#' }
#'@examples
#' wildpop <- simulate_admixture(pop_size = 100,
#'                               number_of_founders = 10,
#'                               total_runtime = 100,
#'                               morgan = 1)
#'
#' ld_results <- calculate_ld(pop = wildpop,
#'                            number_of_markers = 10,
#'                            random_markers = TRUE)
#'
#' plot(ld_results$ld_matrix~ld_results$dist_matrix,
#'      pch = 16,
#'      xlab="Distance between markers",
#'      ylab = "Linkage Disequilibrium")
#' @export
calculate_ld <- function(pop,
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

  ld_matrix   <- matrix(nrow = length(markers), ncol = length(markers), NA)
  rsq_matrix  <- matrix(nrow = length(markers), ncol = length(markers), NA)
  dist_matrix <- matrix(nrow = length(markers), ncol = length(markers), NA)

  for (x in seq_along(markers)) {
    for (y in seq_len(x)) {
      if (x != y) {
        index1 <- c((x - 1) * 2 + 1, (x - 1) * 2 + 2)
        index2 <- c((y - 1) * 2 + 1, (y - 1) * 2 + 2)
        g1 <- all_loci[, index1]
        g2 <- all_loci[, index2]

        ld <- calculate_average_ld(g1, g2)
        ld_matrix[x, y] <- ld$LD
        rsq_matrix[x, y] <- ld$r_sq
        gen_dist <- abs(markers[x] - markers[y])
        dist_matrix[x, y] <- gen_dist
      }
    }
  }

  return(list("ld_matrix" = ld_matrix,
              "rsq_matrix" = rsq_matrix,
              "dist_matrix" = dist_matrix))
}

count_ab <- function(alleles_pos_1, alleles_pos_2, a, b) {
  total_count <- 0
  for (i in 1:2) {
    v1 <- alleles_pos_1[, i] == a
    v2 <- alleles_pos_2[, i] == b
    total_count <- total_count + sum((v1 == v2) & v1 == TRUE)
  }

  return(total_count)
}


#' Calculates the  ld between two alleles
#' @description calculate the average ld between two loci
#' @param alleles_pos_1 alleles at locus 1
#' @param alleles_pos_2 alleles at locus 2
#' @return a list with two entries: LD and r_squared
#' @export
calculate_average_ld <- function(alleles_pos_1, alleles_pos_2) {
  all_alleles <- c(as.vector(alleles_pos_1), as.vector(alleles_pos_2))
  all_alleles <- sort(unique(all_alleles))

  ld <- 0
  r_squared <- 0
  for (i in all_alleles) {
    for (j in all_alleles) {
      p_a_i <- length(which(alleles_pos_1 == i)) / length(alleles_pos_1)
      p_b_j <- length(which(alleles_pos_2 == j)) / length(alleles_pos_2)

      countab <- count_ab(alleles_pos_1, alleles_pos_2, i, j)

      p_a_i_b_j <- countab / (length(alleles_pos_1[, 1]) +
                              length(alleles_pos_1[, 2]))

      if (is.nan(p_a_i_b_j)) p_a_i_b_j <- 0

      d_i_j <- p_a_i_b_j - p_a_i * p_b_j

      d_i_j_max <- 1

      if (d_i_j < 0) {
        d_i_j_max <- min(p_a_i * p_b_j, (1 - p_a_i) * (1 - p_b_j))
      } else {
        d_i_j_max <- min(p_a_i * (1 - p_b_j), (1 - p_a_i) * p_b_j)
      }

      if (d_i_j_max > 0) {
        ld <- ld + p_a_i * p_b_j * abs(d_i_j / d_i_j_max)
        r_i_j <- (d_i_j ^ 2) / (p_a_i * (1 - p_a_i) * p_b_j * (1 - p_b_j))
        r_squared <- r_squared + p_a_i * p_b_j * r_i_j
      }
    }
  }

  return(list("LD" = ld,
              "r_sq" = r_squared)
  )
}
