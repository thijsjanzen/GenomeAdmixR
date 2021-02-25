#' Individual based simulation of the breakdown of contiguous ancestry blocks in
#' two populations linked by migration
#' @description Individual based simulation of the breakdown of contiguous
#' ancestry blocks, with or without selection. Simulations can be started from
#' scratch, or from a predefined input population. Two populations are
#' simulated, connected by migration
#' @param input_population_1 Potential earlier simulated population used as
#' starting point for the simulation. If not provided by the user, the
#' simulation starts from scratch.
#' @param input_population_2 Potential earlier simulated population used as
#' starting point for the simulation. If not provided by the user,
#' the simulation starts from scratch.
#' @param pop_size Vector containing the number of individuals in both
#' populations.
#' @param initial_frequencies A list describing the initial frequency of each
#' ancestor in each population. Each entry in the list contains a vector with
#' the frequencies for all ancestor. The length of the vector indicates the
#' number of unique ancestors. If a vector not summing to 1 is provided, the
#' vector is normalized.
#' @param total_runtime  Number of generations
#' @param morgan Length of the chromosome in Morgan (e.g. the number of
#' crossovers during meiosis)
#' @param num_threads number of threads. Default is 1. Set to -1 to use all
#' available threads
#' @param select_matrix Selection matrix indicating the markers which are under
#' selection. If not provided by the user, the simulation proceeds neutrally.
#' If provided, each row in the matrix should contain five entries:
#' \code{location}{ location of the marker under selection (in Morgan) }
#' \code{fitness of wildtype (aa)} \code{fitness of heterozygote (aA)}
#' \code{fitness of homozygote mutant (AA)} \code{Ancestral type that
#' representes the mutant allele A}
#' @param verbose Verbose output if TRUE. Default value is FALSE
#' @param markers A vector of locations of markers (relative locations in
#' [0, 1]). If a vector is provided, ancestry at these marker positions is
#' tracked for every generation.
#' @param track_junctions Track the average number of junctions over time if
#' TRUE
#' @param multiplicative_selection Default: TRUE. If TRUE, fitness is
#' calculated for multiple markers by multiplying fitness values for each
#' marker. If FALSE, fitness is calculated by adding fitness values for each
#' marker.
#' @param migration_rate  Rate of migration between the two populations.
#' Migration is implemented such that with probability m (migration rate) one
#' of the two parents of a new offspring is from the other population, with
#' probability 1-m both parents are of the focal population.
#' @param stop_at_critical_fst option to stop at a critical FST value
#' , default is FALSE
#' @param critical_fst the critical fst value to stop, if
#' \code{stop_simulation_at_critical_fst} is TRUE
#' @param generations_between_update The number of generations after which the
#' simulation has to check again whether the critical Fst value is exceeded
#' @param sampled_individuals Number of individuals to be sampled at random from
#' the population to estimate Fst
#' @param number_of_markers Number of markers to be used to estimate Fst
#' @param random_markers Are the markers to estimate Fst randomly distributed,
#' or regularly distributed? Default is TRUE.
#' @return A list with: \code{population_1}, \code{population_2} two population
#' objects, and three tibbles with allele frequencies (only contain values of a
#' vector was provided to the argument \code{markers}: \code{frequencies},
#' \code{initial_frequencies} and \code{final_frequencies}. Each tibble contains
#' five columns, \code{time}, \code{location}, \code{ancestor}, \code{frequency}
#' and \code{population}, which indicates the number of generations, the
#' location along the chromosome of the marker, the ancestral allele at that
#' location in that generation, the frequency of that allele and the population
#' in which it was recorded (1 or 2). If a critical fst value was used to
#' terminate the simulation, and object \code{FST} with the final FST estimate
#' is returned as well.
#' @export
simulate_ancestry_migration <- function(input_population_1 = NA,
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
                                         stop_at_critical_fst = FALSE,
                                         critical_fst = 0.1,
                                         generations_between_update = 100,
                                         sampled_individuals = 10,
                                         number_of_markers = 100,
                                         random_markers = TRUE) {

  if (stop_at_critical_fst) {
    message("stopping when FST is: ", critical_fst)
    return(simulate_ancestry_until(input_population_1 = input_population_1,
                                    input_population_2 = input_population_2,
                                    pop_size = pop_size,
                                    initial_frequencies = initial_frequencies,
                                    total_runtime = total_runtime,
                                    morgan = morgan,
                                    num_threads = num_threads,
                                    select_matrix = select_matrix,
                                    markers = markers,
                                    verbose = verbose,
                                    track_junctions = track_junctions,
                                    multiplicative_selection =
                                      multiplicative_selection,
                                    migration_rate = migration_rate,
                                    critical_fst = critical_fst,
                                    generations_between_update =
                                      generations_between_update,
                                    sampled_individuals = sampled_individuals,
                                    number_of_markers = number_of_markers,
                                    random_markers = random_markers))
  }

  if (!is.list(initial_frequencies)) {
    if (length(initial_frequencies) <= 2) {
      stop("need separate frequencies for each population")
    }
  }

  input_population_1 <- check_input_pop(input_population_1)
  input_population_2 <- check_input_pop(input_population_2)

  initial_frequencies <- check_initial_frequencies(initial_frequencies)

  select_matrix <-     check_select_matrix(select_matrix,
                                           markers,
                                           use_data = FALSE)

  if (length(pop_size) == 1) {
    pop_size <- c(pop_size, pop_size)
  }


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

  init_freq_matrix <- matrix(nrow = length(initial_frequencies),
                             ncol = length(initial_frequencies[[1]]))

  for (i in seq_along(initial_frequencies)) {
    for (j in seq_along(initial_frequencies[[i]])) {
      init_freq_matrix[i, j] <- initial_frequencies[[i]][j]
    }
  }

  input_population_1 <- population_to_vector(input_population_1)
  input_population_2 <- population_to_vector(input_population_2)

  selected_pop <- simulate_migration_cpp(input_population_1,
                                         input_population_2,
                                         select_matrix,
                                         pop_size,
                                         init_freq_matrix,
                                         total_runtime,
                                         morgan,
                                         verbose,
                                         track_frequency,
                                         markers,
                                         track_junctions,
                                         multiplicative_selection,
                                         migration_rate,
                                         num_threads)

  selected_popstruct_1 <- create_pop_class(selected_pop$population_1)
  selected_popstruct_2 <- create_pop_class(selected_pop$population_2)

  colnames(selected_pop$initial_frequencies) <- c("time",
                                                  "location",
                                                  "ancestor",
                                                  "frequency",
                                                  "population")

  colnames(selected_pop$final_frequencies) <- c("time",
                                                "location",
                                                "ancestor",
                                                "frequency",
                                                "population")

  initial_freq_tibble <- tibble::as_tibble(selected_pop$initial_frequencies)
  final_freq_tibble   <- tibble::as_tibble(selected_pop$final_frequencies)

  output <- generate_output_list_two_pop(selected_pop,
                                         selected_popstruct_1,
                                         selected_popstruct_2,
                                         initial_freq_tibble,
                                         final_freq_tibble,
                                         track_frequency,
                                         track_junctions)

  class(output) <- "genomadmixr_simulation"
  return(output)
}
