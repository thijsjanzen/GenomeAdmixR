#' create sequence module
#' @description creates a sequence module, which contains all relevant
#' information in order to perform a simulation based on sequence data.
#' @param molecular_data Genomic data used as input, should be of type
#' genomeadmixr_data. Either a single dataset is provided, or a list of
#' multiple genomeadmixr_data objects.
#' @param initial_frequencies A vector describing the initial contribution of
#' each provided input data set to the starting hybrid swarm. By default, equal
#' frequencies are assumed. If a vector not summing to 1 is provided, the vector
#' is normalized.
#' @param morgan Length of the molecular sequence in Morgan (e.g. the number of
#' crossovers during meiosis), alternatively, the recombination rate can be
#' used, see below.
#' @param recombination_rate rate in cM / Mbp, used to map recombination to the
#' markers. If the recombination_rate is not set, the value for Morgan is used,
#' assuming that the markers included span an entire chromosome.
#' @param markers A vector of locations of markers, these markers are
#' tracked for every generation.
#' @param mutation_rate the per base probability of mutation. Default is 0.
#' @param substitution_matrix a 4x4 matrix representing the probability of
#' mutating to another base (where [1/2/3/4] = [a/c/t/g]), conditional on the
#' event of a mutation happening. Default is the JC69 matrix, with equal
#' probabilities for all transitions / transversions.
#' @return sequence module object, used as starting point for
#' simulate_admixture.
#' @export
sequence_module <- function(molecular_data = NA,
                            initial_frequencies = NA,
                            morgan = 1,
                            recombination_rate = NA,
                            markers = NA,
                            mutation_rate = 0,
                            substitution_matrix =
                              matrix(1 / 4, 4, 4)) {

  local_module <- list(input_data = molecular_data,
                       initial_frequencies = initial_frequencies,
                       morgan = morgan,
                       recombination_rate = recombination_rate,
                       markers = markers,
                       mutation_rate = mutation_rate,
                       substitution_matrix = substitution_matrix,
                       type  = "sequence")
  return(local_module)
}
