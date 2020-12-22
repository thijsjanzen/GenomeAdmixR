#' Calculate heterozygosity
#' @description Calculate the average population level heterozygosity
#' @param source_pop Population for which to estimate allele frequencies, or a
#' list of individuals for which to calculate average heterozygosity
#' @param locations A vector indicating the locations (in Morgan) of markers for
#' which to calculate the heterozygosity
#' @param progress_bar Displays a progress_bar if TRUE. Default value is TRUE
#' @return A tibble containing the heterozygosities
#' @export
calculate_heterozygosity <- function(source_pop,
                                     locations,
                                     progress_bar = TRUE) {

  source_pop <- check_input_pop(source_pop)

  pop_for_cpp <- population_to_vector(source_pop)

  heterozygosity_table <- calculate_heterozygosity_cpp(pop_for_cpp,
                                                       locations,
                                                       progress_bar)

  output <- heterozygosity_table
  colnames(output) <- c("location", "heterozygosity")
  output <- tibble::as_tibble(output)

  return(output)
}
