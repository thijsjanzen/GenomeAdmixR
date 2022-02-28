#' Individual based simulation of the breakdown of contiguous ancestry blocks.
#' @description Individual based simulation of the breakdown of contiguous
#' ancestry blocks, with or without selection. Simulations can be started from
#' scratch, or from a predefined input population.
#' @param input_population Potential earlier simulated population used as
#' starting point for the simulation. If not provided by the user, the
#' simulation starts from scratch.
#' @param pop_size The number of individuals in the population. If the number is
#' larger than the number of individuals in the input population (if provided),
#' additional individuals are sampled randomly from the input population to
#' reach the intended size.
#' @param number_of_founders Number of unique ancestors
#' @param initial_frequencies A vector describing the initial frequency of each
#' ancestor. By default, equal frequencies are assumed. If a vector not summing
#' to 1 is provided, the vector is normalized.
#' @param total_runtime  Number of generations
#' @param morgan Length of the chromosome in Morgan (e.g. the number of
#' crossovers during meiosis)
#' @param num_threads number of threads. Default is 1. Set to -1 to use all
#' available threads
#' @param select_matrix Selection matrix indicating the markers which are under
#' selection. If not provided by the user, the simulation proceeds neutrally. If
#' provided, each row in the matrix should contain five entries:
#' \code{location}{ location of the marker under selection (in Morgan) }
#' \code{fitness of wildtype (aa)} \code{fitness of heterozygote (aA)}
#' \code{fitness of homozygote mutant (AA)} \code{Ancestral type that
#' represents the mutant allele A}
#' @param verbose Verbose output if TRUE. Default value is FALSE
#' @param markers A vector of locations of markers (relative locations in
#' [0, 1]). If a vector is provided, ancestry at these marker positions is
#' tracked for every generation.
#' @param track_junctions Track the average number of junctions over time if
#' TRUE
#' @param multiplicative_selection Default: TRUE. If TRUE, fitness is calculated
#' for multiple markers by multiplying fitness values for each marker. If FALSE,
#' fitness is calculated by adding fitness values for each marker.
#' @return A list with: \code{population} a population object, and three tibbles
#' with allele frequencies (only contain values of a vector was provided to the
#' argument \code{markers}: \code{frequencies} , \code{initial_frequencies} and
#' \code{final_frequencies}. Each tibble contains four columns, \code{time},
#' \code{location}, \code{ancestor} and \code{frequency}, which indicates the
#' number of generations, the location along the chromosome of the marker, the
#' ancestral allele at that location in that generation, and finally, the
#' frequency of that allele.
simulate_ancestry <- function(input_population = NA,
                               pop_size = NA,
                               number_of_founders = 2,
                               initial_frequencies = NA,
                               total_runtime = 100,
                               morgan = 1,
                               num_threads = 1,
                               select_matrix = NA,
                               markers = NA,
                               verbose = FALSE,
                               track_junctions = FALSE,
                               multiplicative_selection = TRUE) {



  input_population <- check_input_pop(input_population)

  if (is.na(pop_size)) {
    if (inherits(input_population, "population")) {
      pop_size <- length(input_population)
    } else {
      stop("pop_size is undefined, need an input population")
    }
  }

  select_matrix <- check_select_matrix(select_matrix,
                                       markers,
                                       use_data = FALSE)

  if (length(markers) == 1) {
    if (is.na(markers))  {
      markers <- c(-1, -1)
      track_frequency <- FALSE
    } else {
      track_frequency <- TRUE
    }
  } else {
    track_frequency <- TRUE
  }

  if (inherits(input_population, "population")) {
    input_population <- population_to_vector(input_population)
  }

  selected_pop <- simulate_cpp(input_population,
                               select_matrix,
                               pop_size,
                               number_of_founders,
                               initial_frequencies,
                               total_runtime,
                               morgan,
                               verbose,
                               track_frequency,
                               markers,
                               track_junctions,
                               multiplicative_selection,
                               num_threads)

  selected_popstruct <- create_pop_class(selected_pop$population)

  colnames(selected_pop$initial_frequencies) <- c("time",
                                                  "location",
                                                  "ancestor",
                                                  "frequency")
  colnames(selected_pop$final_frequencies) <-
     colnames(selected_pop$initial_frequencies)


  initial_freq_tibble <- tibble::as_tibble(selected_pop$initial_frequencies)
  final_freq_tibble   <- tibble::as_tibble(selected_pop$final_frequencies)


  output <- generate_output_list_one_pop(selected_popstruct,
                                         selected_pop,
                                         initial_freq_tibble,
                                         final_freq_tibble,
                                         track_frequency,
                                         track_junctions)

  class(output) <- "genomadmixr_simulation"

  return(output)
}
