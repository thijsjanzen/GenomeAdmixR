#' Calculate FST
#' @description The FST value between two populations is calculated, given a
#' number of markers. Markers are superimposed upon the (known) ancestry along
#' the chromosome for all sampled individuals. Markers can be chosen to be
#' regularly spaced, or randomly distributed.
#' @param pop1 Population object
#' @param pop2 Population object
#' @param sampled_individuals Number of individuals to base the FST upon.
#' Individuals are randomly drawn from each population, a lower number speeds
#' up calculations.
#' @param number_of_markers Number of markers along the chromosome used to
#' calculate FST metrics.
#' @param random_markers  If TRUE, markers are randomly spaced along the
#' chromosome, if FALSE, markers are equidistantly spaced along the chromosome.
#' @details Uses the function \code{wc} from the package \code{hierfstat} to
#' calculate the FST. The function \code{wc} computes the Weir and Cockerham
#' F statistic.
#' @return FST value
#' @examples
#' two_populations <- simulate_admixture_migration(pop_size = c(100, 100),
#'                                                 morgan = 1,
#'                                                 migration_rate = 0.01)
#'
#' FST <- calculate_fst(pop1 = two_populations$population_1,
#'                      pop2 = two_populations$population_2,
#'                      sampled_individuals = 10,
#'                      number_of_markers = 100,
#'                      random_markers = TRUE)
#' @export
calculate_fst <- function(pop1,
                          pop2,
                          sampled_individuals,
                          number_of_markers = 100,
                          random_markers = FALSE) {

  pop1 <- check_input_pop(pop1)

  pop2 <- check_input_pop(pop2)

  number_of_markers <- round(number_of_markers)

  all_loci <- create_loci_matrix(
                pop1[sample(seq_along(pop1), sampled_individuals)],
                pop2[sample(seq_along(pop2), sampled_individuals)],
                number_of_markers,
                random_markers)

  hierf_wc <- hierfstat::wc(as.data.frame(all_loci))

  return(hierf_wc$FST)
}


#' @keywords internal
create_loci_matrix <- function(pop1,
                               pop2,
                               number_of_markers,
                               random_markers) {

  all_loci <- matrix(nrow = length(pop1) + length(pop2),
                     ncol = 1 + number_of_markers, 0)
  all_loci[, 1] <- c(rep(1, length(pop1)), rep(2, length(pop1)))
  colnames(all_loci) <- c("population", 1:number_of_markers)

  markers <- seq(1e-9, 1 - (1e-9), length.out = number_of_markers)
  if (random_markers) {
    markers <- create_random_markers(number_of_markers)
  }

  for (x in seq_along(markers)) {
    focal_marker <- markers[x]
    for (i in seq_along(pop1)) {
      allele_1 <- 10 + findtype(pop1[[i]]$chromosome1, focal_marker)
      allele_2 <- 10 + findtype(pop1[[i]]$chromosome2, focal_marker)
      final_allele <- paste0(allele_1, allele_2)
      all_loci[i, x + 1] <- as.numeric(final_allele)
    }

    for (i in seq_along(pop2)) {
      allele_1 <- 10 + findtype(pop2[[i]]$chromosome1, focal_marker)
      allele_2 <- 10 + findtype(pop2[[i]]$chromosome2, focal_marker)
      final_allele <- paste0(allele_1, allele_2)
      all_loci[length(pop1) + i, x + 1] <- as.numeric(final_allele)
    }
  }

  return(all_loci)
}
