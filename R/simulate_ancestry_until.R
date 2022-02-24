#' @keywords internal
simulate_ancestry_until <- function(input_population_1 = NA,
                                     input_population_2 = NA,
                                     pop_size = c(100, 100),
                                     initial_frequencies = list(c(1.0, 0),
                                                                c(0, 1.0)),
                                     total_runtime = 100,
                                     morgan = 1,
                                     num_threads = 1,
                                     select_matrix = NA,
                                     markers = NA,
                                     verbose = FALSE,
                                     track_junctions = FALSE,
                                     multiplicative_selection = TRUE,
                                     migration_rate = 0.0,
                                     generations_between_update = 100,
                                     critical_fst = 0.1,
                                     sampled_individuals = 10,
                                     number_of_markers = 100,
                                     random_markers = TRUE) {

  pops <- simulate_ancestry_migration(
    input_population_1 = input_population_1,
    input_population_2 = input_population_2,
    pop_size = pop_size,
    initial_frequencies = initial_frequencies,
    total_runtime = generations_between_update,
    num_threads = num_threads,
    morgan = morgan,
    select_matrix = select_matrix,
    markers = markers,
    verbose = verbose,
    track_junctions = track_junctions,
    multiplicative_selection = multiplicative_selection,
    migration_rate = migration_rate)

  fst <- calculate_fst(pop1 = pops$population_1,
                       pop2 = pops$population_2,
                       sampled_individuals = sampled_individuals,
                       number_of_markers = number_of_markers,
                       random_markers = random_markers)

  if (verbose) message("Number of Generations\tFST")
  if (verbose) message(generations_between_update, "\t", fst)

  total_generations <- generations_between_update
  while (fst < critical_fst && total_generations < total_runtime) {

    pops <- simulate_ancestry_migration(
      input_population_1 = pops$population_1,
      input_population_2 = pops$population_2,
      pop_size = pop_size,
      initial_frequencies = initial_frequencies,
      total_runtime = generations_between_update,
      morgan = morgan,
      num_threads = num_threads,
      select_matrix = select_matrix,
      markers = markers,
      verbose = verbose,
      track_junctions = track_junctions,
      multiplicative_selection = multiplicative_selection,
      migration_rate = migration_rate)

    fst <- calculate_fst(pops$population_1, pops$population_2,
                         sampled_individuals = sampled_individuals,
                         number_of_markers = number_of_markers,
                         random_markers = random_markers)

    total_generations <- total_generations + generations_between_update
    if (verbose) message(total_generations, "\t", fst)
  }

  return(list("population_1" = pops$population_1,
              "population_2" = pops$population_2,
              "Number_of_generations" = total_generations,
              "FST" = fst))
}
