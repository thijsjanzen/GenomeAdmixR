#' Individual based simulation of the breakdown of contiguous ancestry blocks.
#' @description Individual based simulation of the breakdown of contiguous
#' ancestry blocks, with or without selection. Simulations can be started from
#' scratch, or from a predefined input population.
#' @param input_data Genomic data used as input, should be of type
#' genomeadmixr_data. Either a single dataset is provided, or a list of
#' multiple genomeadmixr_data objects.
#' @param pop_size Vector containing the number of individuals in both
#' populations.
#' @param initial_frequencies A vector describing the initial contribution of
#' each provided input data set to the starting hybrid swarm. By default, equal
#' frequencies are assumed. If a vector not summing to 1 is provided, the vector
#' is normalized.
#' @param total_runtime  Number of generations
#' @param morgan Length of the chromosome in Morgan (e.g. the number of
#' crossovers during meiosis)
#' @param recombination_rate rate in cM / Mbp, used to map recombination to the
#' markers. If the recombination_rate is not set, the value for Morgan is used,
#' assuming that the markers included span an entire chromosome.
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
#' @param markers A vector of locations of markers, these markers are
#' tracked for every generation.
#' @param multiplicative_selection Default: TRUE. If TRUE, fitness is calculated
#' for multiple markers by multiplying fitness values for each marker. If FALSE,
#' fitness is calculated by adding fitness values for each marker.
#' @param mutation_rate the per base probability of mutation. Default is 0.
#' @param substitution_matrix a 4x4 matrix representing the probability of
#' mutating to another base (where [1/2/3/4] = [a/c/t/g]), conditional on the
#' event of a mutation happening. Default is the JC69 matrix, with equal
#' probabilities for all transitions / transversions.
#' @return A list with: \code{population} a population object, and three tibbles
#' with allele frequencies (only contain values of a vector was provided to the
#' argument \code{markers}: \code{frequencies} , \code{initial_frequencies} and
#' \code{final_frequencies}. Each tibble contains four columns, \code{time},
#' \code{location}, \code{ancestor} and \code{frequency}, which indicates the
#' number of generations, the location along the chromosome of the marker, the
#' ancestral allele at that location in that generation, and finally, the
#' frequency of that allele.
simulate_sequence <- function(input_data = NA,
                              pop_size = NA,
                              initial_frequencies = NA,
                              total_runtime = 100,
                              morgan = 1,
                              recombination_rate = NA,
                              num_threads = 1,
                              select_matrix = NA,
                              markers = NA,
                              verbose = FALSE,
                              multiplicative_selection = TRUE,
                              mutation_rate = 0,
                              substitution_matrix = matrix(1 / 4, 4, 4)) {

  input_data <- verify_genomeadmixr_data(input_data,
                                         markers, verbose = verbose)

  if (!inherits(input_data, "genomeadmixr_data")) {
    if (is.list(input_data)) {
      if (length(input_data) > 1) {
        if (verbose) message("found multiple input populations")
        input_data <- combine_input_data(input_data,
                                         frequencies = initial_frequencies,
                                         pop_size = pop_size)
      }
    }
  }

  if (is.na(pop_size)) {
    if (length(input_data) == 2) {
      pop_size <- length(input_data$genomes[, 1]) / 2
    } else {
      stop("pop_size is undefined, need an input population")
    }
  }

  select_matrix <- check_select_matrix(select_matrix,
                                       markers,
                                       use_data = TRUE)

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

  if (track_frequency) {
    markers <- check_markers(markers, input_data$markers)
  }

  if (mutation_rate > 0) {
    substitution_matrix <- verify_substitution_matrix(substitution_matrix)
  }

  recombination_map <- c(-1, 0)
  if (!is.na(recombination_rate)) {
    recombination_map <- create_recombination_map(input_data$markers,
                                                  recombination_rate)
  }

  selected_pop <- simulate_emp_cpp(input_data$genomes,
                                   input_data$markers,
                                   select_matrix,
                                   pop_size,
                                   total_runtime,
                                   morgan,
                                   verbose,
                                   track_frequency,
                                   markers,
                                   multiplicative_selection,
                                   mutation_rate,
                                   substitution_matrix,
                                   num_threads,
                                   recombination_map)

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
                                         track_frequency)

  output$frequencies <- convert_to_dna(output$frequencies)
  output$initial_frequency <- convert_to_dna(output$initial_frequency)
  output$final_frequency <- convert_to_dna(output$final_frequency)

  class(output) <- "genomadmixr_simulation"

  return(output)
}

#' @keywords internal
convert_to_dna <- function(frequency_tibble) {
  alleles <- frequency_tibble$ancestor
  alleles[alleles == 1] <- "a"
  alleles[alleles == 2] <- "c"
  alleles[alleles == 3] <- "t"
  alleles[alleles == 4] <- "g"
  alleles[alleles == 0] <- "-"

  frequency_tibble$ancestor <- alleles
  return(frequency_tibble)
}
