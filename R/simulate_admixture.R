#' Individual based simulation of the breakdown of contiguous ancestry blocks.
#' @description Individual based simulation of the breakdown of contiguous
#' ancestry blocks, with or without selection. Simulations can be started from
#' scratch, or from a predefined input population.
#' @param module Chosen module to simulate, either created with
#' \code{module_ancestry} or
#' \code{module_sequence}.
#' @param pop_size The number of individuals in the population. If the number is
#' larger than the number of individuals in the input population (if provided),
#' additional individuals are sampled randomly from the input population to
#' reach the intended size.
#' @param total_runtime  Number of generations
#' @param select_matrix Selection matrix indicating the markers which are under
#' selection. If not provided by the user, the simulation proceeds neutrally. If
#' provided, each row in the matrix should contain five entries:
#' \code{location}{ location of the marker under selection (in Morgan) }
#' \code{fitness of wildtype (aa)} \code{fitness of heterozygote (aA)}
#' \code{fitness of homozygote mutant (AA)} \code{Ancestral type that
#' represents the mutant allele A}
#' @param multiplicative_selection Default: TRUE. If TRUE, fitness is calculated
#' for multiple markers by multiplying fitness values for each marker. If FALSE,
#' fitness is calculated by adding fitness values for each marker.
#' @param verbose Verbose output if TRUE. Default value is FALSE
#' @param num_threads number of threads. Default is 1. Set to -1 to use all
#' available threads
#' @return A list with: \code{population} a population object, and three tibbles
#' with allele frequencies (only contain values of a vector was provided to the
#' argument \code{markers}: \code{frequencies} , \code{initial_frequencies} and
#' \code{final_frequencies}. Each tibble contains four columns, \code{time},
#' \code{location}, \code{ancestor} and \code{frequency}, which indicates the
#' number of generations, the location along the chromosome of the marker, the
#' ancestral allele at that location in that generation, and finally, the
#' frequency of that allele.
#' @export
simulate_admixture <- function(module = ancestry_module(),
                               pop_size = 100,
                               total_runtime = 100,
                               select_matrix = NA,
                               multiplicative_selection = TRUE,
                               verbose = FALSE,
                               num_threads = 1) {

  if (is.na(module$migration$migration_rate)) {
    if (module$type == "ancestry") {
      result <- simulate_ancestry(input_population = module$input_data,
                                  pop_size = pop_size,
                                  number_of_founders = module$number_of_founders,
                                  initial_frequencies =
                                    module$initial_frequencies,
                                  total_runtime = total_runtime,
                                  morgan = module$morgan,
                                  num_threads = num_threads,
                                  select_matrix = select_matrix,
                                  markers = module$markers,
                                  verbose = verbose,
                                  track_junctions = module$track_junctions,
                                  multiplicative_selection =
                                    multiplicative_selection)
      return(result)
    }

    if (module$type == "sequence") {
      result <- simulate_sequence(input_data = module$input_data,
                                  pop_size = pop_size,
                                  initial_frequencies =
                                    module$initial_frequencies,
                                  total_runtime = total_runtime,
                                  morgan = module$morgan,
                                  recombination_rate = module$recombination_rate,
                                  num_threads = num_threads,
                                  select_matrix = select_matrix,
                                  markers = module$markers,
                                  verbose = verbose,
                                  multiplicative_selection =
                                    multiplicative_selection,
                                  mutation_rate = module$mutation_rate,
                                  substitution_matrix =
                                    module$substitution_matrix)
      return(result)
    }
  } else {
    if (verbose)
      message("found positive migration rate, assuming two connected populations")
    if (module$type == "ancestry") {
      result <- simulate_ancestry_migration(input_population_1 =
                                              module$input_data$population_1,
                                            input_population_2 =
                                              module$input_data$population_2,
                                            pop_size =
                                              module$migration$population_size,
                                            initial_frequencies =
                                              module$migration$initial_frequencies,
                                            total_runtime = total_runtime,
                                            morgan = module$morgan,
                                            num_threads = num_threads,
                                            select_matrix = select_matrix,
                                            markers = module$markers,
                                            verbose = verbose,
                                            track_junctions =
                                              module$track_junctions,
                                            multiplicative_selection =
                                              multiplicative_selection,
                                            migration_rate =
                                              module$migration$migration_rate,
                                            stop_at_critical_fst =
                                              module$migration$stop_at_critical_fst,
                                            critical_fst =
                                              module$migration$critical_fst,
                                            generations_between_update =
                                              module$migration$generations_between_update,
                                            sampled_individuals =
                                              module$migration$sampled_individuals,
                                            number_of_markers =
                                              module$migration$number_of_markers,
                                            random_markers =
                                              module$migration$random_markers)

      return(result)
    }
    if (module$type == "sequence") {
      if (verbose)
        message("found positive migration rate, assuming two connected populations")
      result <- simulate_sequence_migration(input_data_population_1 =
                                              module$input_data$population_1,
                                            input_data_population_2 =
                                              module$input_data$population_2,
                                            pop_size =
                                              module$migration$population_size,
                                            total_runtime = total_runtime,
                                            morgan = module$morgan,
                                            recombination_rate =
                                              module$recombination_rate,
                                            num_threads = num_threads,
                                            select_matrix = select_matrix,
                                            markers = module$markers,
                                            verbose = verbose,
                                            multiplicative_selection =
                                              multiplicative_selection,
                                            migration_rate =
                                              module$migration$migration_rate,
                                            stop_at_critical_fst =
                                              module$migration$stop_at_critical_fst,
                                            critical_fst = module$migration$critical_fst,
                                            generations_between_update =
                                              module$migration$generations_between_update,
                                            sampled_individuals =
                                              module$migration$sampled_individuals,
                                            number_of_markers =
                                              module$migration$number_of_markers,
                                            random_markers =
                                              module$migration$random_markers,
                                            mutation_rate =
                                              module$mutation_rate,
                                            substitution_matrix =
                                              module$substitution_matrix)

      return(result)
    }
  }
  stop("wrong module used, please pick a module with
       ancestry_module() or sequence_module()")
}
