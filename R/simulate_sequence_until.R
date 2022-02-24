#' @keywords internal
simulate_sequence_until <- function(input_data_population_1 = NA,
                                    input_data_population_2 = NA,
                                    pop_size = c(100, 100),
                                    total_runtime = 100,
                                    morgan = 1,
                                    recombination_rate = NA,
                                    num_threads = 1,
                                    select_matrix = NA,
                                    markers = NA,
                                    verbose = FALSE,
                                    multiplicative_selection = TRUE,
                                    migration_rate = 0.0,
                                    generations_between_update = 100,
                                    critical_fst = 0.1,
                                    sampled_individuals = 10,
                                    number_of_markers = 100,
                                    random_markers = TRUE,
                                    mutation_rate = 0,
                                    substitution_matrix =
                                      matrix(1 / 4, nrow = 4, ncol = 4)) {

  pops <- simulate_sequence_migration(
    input_data_population_1 = input_data_population_1,
    input_data_population_2 = input_data_population_2,
    pop_size = pop_size,
    total_runtime = generations_between_update,
    num_threads = num_threads,
    morgan = morgan,
    recombination_rate = recombination_rate,
    select_matrix = select_matrix,
    markers = markers,
    verbose = verbose,
    multiplicative_selection = multiplicative_selection,
    migration_rate = migration_rate,
    mutation_rate = mutation_rate,
    substitution_matrix = substitution_matrix)

  fst <- calculate_fst(pops$population_1, pops$population_2,
                       sampled_individuals = sampled_individuals,
                       number_of_markers = number_of_markers,
                       random_markers = random_markers)

  if (verbose) message("Number of Generations\tFST\n")
  if (verbose) message(generations_between_update, "\t", fst, "\n")

  total_generations <- generations_between_update
  while (fst < critical_fst && total_generations < total_runtime) {

    pop1_for_data_cpp <-
      simulation_data_to_genomeadmixr_data(pops$population_1,
                                           input_data_population_1$markers)
    pop2_for_data_cpp <-
      simulation_data_to_genomeadmixr_data(pops$population_2,
                                           input_data_population_1$markers)

    pops <- simulate_sequence_migration(
      input_data_population_1 = pop1_for_data_cpp,
      input_data_population_2 = pop2_for_data_cpp,
      pop_size = pop_size,
      total_runtime = generations_between_update,
      morgan = morgan,
      recombination_rate = recombination_rate,
      num_threads = num_threads,
      select_matrix = select_matrix,
      markers = markers,
      verbose = verbose,
      multiplicative_selection = multiplicative_selection,
      migration_rate = migration_rate,
      mutation_rate = mutation_rate,
      substitution_matrix = substitution_matrix)

    fst <- calculate_fst(pops$population_1, pops$population_2,
                         sampled_individuals = sampled_individuals,
                         number_of_markers = number_of_markers,
                         random_markers = random_markers)

    total_generations <- total_generations + generations_between_update
    message(total_generations, "\t", fst, "\n")
  }
  return(list("population_1" = pops$population_1,
              "population_2" = pops$population_2,
              "Number_of_generations" = total_generations,
              "FST" = fst))
}
