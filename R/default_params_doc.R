#' Default parameter documentation
#'
#' This function's purpose is to list all parameter documentation to be
#' inherited by the relevant functions.
#'
#' @param module Chosen module to simulate, either created with
#' \code{module_ancestry} or
#' \code{module_sequence}.
#' @param pop_size The number of individuals in the population. If the number is
#' larger than the number of individuals in the input population (if provided),
#' additional individuals are sampled randomly from the input population to
#' reach the intended size.
#' @param total_runtime  Number of generations
#' @param migration settings associated with migration, should be created with
#' \code{\link{migration_settings}}
#' @param select_matrix Selection matrix indicating the markers which are under
#' selection. If not provided by the user, the simulation proceeds neutrally. If
#' provided, each row in the matrix should contain five entries:
#' \itemize{
#' \item{location of the marker under selection (in Morgan) }
#' \item{fitness of wildtype (aa)}
#' \item{fitness of heterozygote (aA)}
#' \item{fitness of homozygote mutant (AA)}
#' \item{Ancestral type that represents the mutant allele A}
#' }
#' @param multiplicative_selection Default: TRUE. If TRUE, fitness is calculated
#' for multiple markers by multiplying fitness values for each marker. If FALSE,
#' fitness is calculated by adding fitness values for each marker.
#' @param verbose Verbose output if TRUE. Default value is FALSE
#' @param num_threads number of threads. Default is 1. Set to -1 to use all
#' available threads
#' @param morgan Length of the chromosome in Morgan (e.g. the number of
#' crossovers during meiosis)
#' @param markers A vector of locations of markers (relative locations in
#' [0, 1]). If a vector is provided, ancestry at these marker positions is
#' tracked for every generation.
#' @param track_junctions Track the average number of junctions over time if
#' TRUE
#' @return Nothing
#' @keywords internal
#' @export
default_params_doc <- function(module,
                               total_runtime,
                               migration,
                               select_matrix,
                               pop_size,
                               multiplicative_selection,
                               verbose,
                               num_threads,
                               morgan,
                               markers,
                               track_junctions) {
  # Nothing
}