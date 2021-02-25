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
#' @param migration settings associated with migration, should be created with
#' \code{\link{migration_settings}}
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
#' @examples
#' # local ancestry simulation
#' two_populations <- simulate_admixture(
#'                          module = ancestry_module(number_of_founders = 3,
#'                                                   morgan = 0.8),
#'                          migration = migration_settings(
#'                                          migration_rate = 0.01,
#'                                          population_size = c(100, 100)),
#'                          total_runtime = 10)
#'  # sequence simulation
#'  data(dgrp2.3R.5k.data)
#'
#' sequence_population <-
#'       simulate_admixture(
#'                   module = sequence_module(molecular_data = dgrp2.3R.5k.data,
#'                            recombination_rate = 0.2,
#'                            mutation_rate = 1e-5),
#'                   pop_size = 1000,
#'                   total_runtime = 10)
#' @export
simulate_admixture <- function(module = ancestry_module(),
                               pop_size = 100,
                               total_runtime = 100,
                               migration = migration_settings(),
                               select_matrix = NA,
                               multiplicative_selection = TRUE,
                               verbose = FALSE,
                               num_threads = 1) {

  if (is.na(migration$migration_rate)) {
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

      if (!methods::is(module$input_population, "genomadmixr_simulation")) {
        if (is.list(module$input_population)) {
          input_population2 <- list()
          input_population2$population_1 <- module$input_population[[1]]
          input_population2$population_2 <- module$input_population[[2]]
          module$input_population <- input_population2
        } else {
          if (length(module$input_population) == 1) {
            if (is.na(module$input_population)) {
              module$input_population <- list(population_1 = c(-1e6, -1e6),
                                              population_2 = c(-1e6, -1e6))
            }
          }
        }
      }

      if (verbose) {
        message("starting simulate_ancestry_migration")
        Sys.sleep(1)
      }

      result <- simulate_ancestry_migration(input_population_1 =
                                              module$input_population$population_1,
                                            input_population_2 =
                                              module$input_population$population_2,
                                            pop_size =
                                              migration$population_size,
                                            initial_frequencies =
                                              migration$initial_frequencies,
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
                                              migration$migration_rate,
                                            stop_at_critical_fst =
                                              migration$stop_at_critical_fst,
                                            critical_fst =
                                              migration$critical_fst,
                                            generations_between_update =
                                              migration$generations_between_update,
                                            sampled_individuals =
                                              migration$sampled_individuals,
                                            number_of_markers =
                                              migration$number_of_markers,
                                            random_markers =
                                              migration$random_markers)

      return(result)
    }
    if (module$type == "sequence") {

      if (!methods::is(module$input_data, "genomadmixr_simulation")) {
        if (is.list(module$input_data)) {
          input_population2 <- list()
          input_population2$population_1 <- module$input_data[[1]]
          input_population2$population_2 <- module$input_data[[2]]
          module$input_data <- input_population2
        } else {
          if (length(module$input_data) == 1) {
            if (is.na(module$input_data)) {
              module$input_data <- list(population_1 = NA,
                                        population_2 = NA)
            }
          }
        }
      }

      result <- simulate_sequence_migration(input_data_population_1 =
                                              module$input_data$population_1,
                                            input_data_population_2 =
                                              module$input_data$population_2,
                                            pop_size =
                                              migration$population_size,
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
                                              migration$migration_rate,
                                            stop_at_critical_fst =
                                              migration$stop_at_critical_fst,
                                            critical_fst =
                                              migration$critical_fst,
                                            generations_between_update =
                                              migration$generations_between_update,
                                            sampled_individuals =
                                              migration$sampled_individuals,
                                            number_of_markers =
                                              migration$number_of_markers,
                                            random_markers =
                                              migration$random_markers,
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
