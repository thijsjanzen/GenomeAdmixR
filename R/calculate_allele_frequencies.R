#' Calculate allele frequencies
#' @description Calculate for a number of regularly spaced markers the relative
#' frequency of each ancestor in the population.
#' @param source_pop Population for which to estimate allele frequencies
#' @param locations A vector indicating the locations (in Morgan) where
#' to calculate the allele frequencies.
#' @param progress_bar Displays a progress_bar if TRUE. Default value is TRUE
#' @details Markers are equidistantly spaced, with a distance of
#' \code{step_size} in between them.
#' @return A tibble containing the allele frequencies
#' @examples
#' number_founders = 20
#' wildpop =  simulate_admixture(pop_size = 1000,
#'                               number_of_founders = number_founders,
#'                               total_runtime = 10,
#'                               morgan = 1)
#'
#' freq_output <- calculate_allele_frequencies(wildpop,
#'                                             progress_bar = TRUE)
#'
#' require(ggplot2)
#' ggplot(freq_output, aes(x=location, y = frequency,
#'                         col = as.factor(ancestor))) +
#'   geom_line()
#' @export
calculate_allele_frequencies <- function(source_pop,
                                         locations = seq(0, 1,
                                                         length.out = 100),
                                         progress_bar = TRUE) {

  source_pop <- check_input_pop(source_pop)

  pop_for_cpp <- population_to_vector(source_pop)

  frequency_table <- calculate_allele_spectrum_cpp(pop_for_cpp,
                                                   locations,
                                                   progress_bar)

  output <- frequency_table
  colnames(output) <- c("time", "location", "ancestor", "frequency")
  output <- tibble::as_tibble(output)

  return(output)
}
