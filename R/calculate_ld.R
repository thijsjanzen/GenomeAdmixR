#' Calculate linkage disequilibrium statistics
#' This function calculates two matrices, once containing all pairwise
#' linkage disequilibrium (ld) values, and one matrix containing all pairwise r
#' statistics
#' @param pop focal population
#' @param sampled_individuals Number of individuals randomly sampled to
#' calculate the LD matrices
#' @param markers vector of markers. If a single number is used, that number of
#' markers is randomly placed along the genome.
#' @param verbose display verbose output, default is FALSE.
#' @return An object containing two items:
#' \item{ld_matrix}{
#'   Pairwise ld statistics for all markers
#' }
#' \item{rsq_matrix}{
#'   Pairwise rsq statistics for all markers
#' }
#'@examples
#' wildpop =  simulate_admixture(
#'    module = ancestry_module(number_of_founders = 10, morgan = 1),
#'    pop_size = 1000,
#'    total_runtime = 10)
#'
#' ld_results <- calculate_ld(pop = wildpop,
#'                            markers = 10)
#'
#' plot(ld_results$ld_matrix~ld_results$dist_matrix,
#'      pch = 16,
#'      xlab="Distance between markers",
#'      ylab = "Linkage Disequilibrium")
#' @export
calculate_ld <- function(pop,
                         sampled_individuals = 10,
                         markers = NA,
                         verbose = FALSE) {

  pop <- check_input_pop(pop)

  if (length(markers) == 1) {
    if (is.na(markers)) {
      if (inherits(pop, "genomadmixr_data")) {
        markers <- pop$markers
      }
    } else {
      markers <- create_random_markers(markers[1])
    }
  }

  popx <- pop[1:sampled_individuals]
  class(popx) <- "population"
  pop_for_cpp <- population_to_vector(popx)
  if (verbose) message("starting creation of allele matrix")
  all_loci <- simulation_data_to_genomeadmixr_data_cpp(pop_for_cpp,
                                                        markers)

  all_loci[all_loci < 0] <- NA

  ld_matrix   <- matrix(nrow = length(markers), ncol = length(markers), NA)
  rsq_matrix  <- matrix(nrow = length(markers), ncol = length(markers), NA)
  dist_matrix <- matrix(nrow = length(markers), ncol = length(markers), NA)

  if (verbose) message("done creation of allele matrix")

  if (verbose) pb <- utils::txtProgressBar(min = 0,
                                           max = length(markers),
                                           style = 3)
  for (x in seq_along(markers)) {
    for (y in seq_len(x)) {
      if (x != y) {
        g1 <- all_loci[, x]
        g2 <- all_loci[, y]

        ld <- calculate_average_ld(g1, g2)
        ld_matrix[x, y] <- ld$LD
        rsq_matrix[x, y] <- ld$r_sq
        gen_dist <- abs(markers[x] - markers[y])
        dist_matrix[x, y] <- gen_dist
      }
    }
    if (verbose) utils::setTxtProgressBar(pb, x)
  }

  return(list("ld_matrix" = ld_matrix,
              "rsq_matrix" = rsq_matrix,
              "dist_matrix" = dist_matrix))
}

count_ab <- function(alleles_pos_1, alleles_pos_2, a, b) {
  v1 <- alleles_pos_1 == a
  v2 <- alleles_pos_2 == b
  total_count <- sum((v1 == v2) & v1 == TRUE)
  return(total_count)
}


#' Calculates the  ld between two alleles
#' @description calculate the average ld between two loci
#' @param alleles_pos_1 alleles at locus 1
#' @param alleles_pos_2 alleles at locus 2
#' @return a list with two entries: LD and r_squared
#' @export
calculate_average_ld <- function(alleles_pos_1,
                                 alleles_pos_2) {
  all_alleles <- c(as.vector(alleles_pos_1), as.vector(alleles_pos_2))

  if (sum(is.na(all_alleles)) > 0) {
    return(list("LD" = NA,
                "r_sq" = NA))
  }

  all_alleles <- sort(unique(all_alleles))

  ld <- 0
  r_squared <- 0
  for (i in all_alleles) {
    for (j in all_alleles) {
      p_a_i <- length(which(alleles_pos_1 == i)) / length(alleles_pos_1)
      p_b_j <- length(which(alleles_pos_2 == j)) / length(alleles_pos_2)

      countab <- count_ab(alleles_pos_1, alleles_pos_2, i, j)

      p_a_i_b_j <- countab / (length(alleles_pos_1))

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
