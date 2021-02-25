#' Creates a module to start simulations tracking local ancestry
#' @description Module to perform simulations based on local ancestry
#' @param input_population Potential earlier simulated population used as
#' starting point for the simulation. If not provided by the user, the
#' simulation starts from scratch.
#' @param number_of_founders Number of unique ancestors / ancestries to be
#' tracked in the simulation
#' @param initial_frequencies A vector describing the initial frequency of each
#' ancestor / ancestry. By default, equal frequencies are assumed. If a vector
#' not summing to 1 is provided, the vector is normalized.
#' @param morgan Length of the genomic stretch simulated, expressed in Morgan
#' (e.g. the number of crossovers during meiosis)
#' @param markers A vector of locations of markers, with the location in Morgan.
#' Ancestry at these marker positions is tracked for every generation.
#' @param track_junctions Tracks the average number of junctions over time if
#' TRUE
#' @export
ancestry_module <- function(input_population = NA,
                            number_of_founders = 2,
                            initial_frequencies = NA,
                            morgan = 1,
                            markers = NA,
                            track_junctions = FALSE) {

  if (sum(is.na(initial_frequencies))) {
    initial_frequencies <- rep(1.0 / number_of_founders,
                               times = number_of_founders)
  }

  if (sum(initial_frequencies) != 1) {
    initial_frequencies <- initial_frequencies / sum(initial_frequencies)
  }

  local_module <- list(input_population = input_population,
                       number_of_founders = number_of_founders,
                       initial_frequencies = initial_frequencies,
                       morgan = morgan,
                       markers = markers,
                       track_junctions = track_junctions,
                       type = "ancestry")

  return(local_module)
}