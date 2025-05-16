#' Individual based simulation of the breakdown of contiguous ancestry blocks in
#' two populations linked by migration
#' @description Individual based simulation of the breakdown of contiguous
#' ancestry blocks, with or without selection. Simulations can be started from
#' scratch, or from a predefined input population. Two populations are
#' simulated, connected by migration
#' @param input_data_population_1 Genomic data used as input, should be created
#' by the function \code{create_input_data} or by the function
#' \code{combine_input_data}
#' @param input_data_population_2 Genomic data used as input, should be created
#' by thefunction \code{create_input_data} or by the function
#' \code{combine_input_data}
#' @param pop_size Vector containing the number of individuals in both
#' populations.
#' @inheritParams default_params_doc
#' @param recombination_rate rate in cM / Mbp, used to map recombination to the
#' markers. If the recombination_rate is not set, the value for morgan is used,
#' assuming that the markers included span an entire chromosome.
#' @param num_threads number of threads. Default is 1. Set to -1 to use all
#' available threads
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
#' @param mutation_rate the per base probability of mutation. Default is 0.
#' @param substitution_matrix a 4x4 matrix representing the probability of
#' mutating to another base (where [1/2/3/4] = [a/c/t/g]), conditional on the
#' event of a mutation happening. Default is the JC69 matrix, with equal
#' probabilities for all transitions / transversions.
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
simulate_sequence_migration <- function(input_data_population_1 = NA, # nolint
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
                                        stop_at_critical_fst = FALSE,
                                        critical_fst = NA,
                                        generations_between_update = 100,
                                        sampled_individuals = 10,
                                        number_of_markers = 100,
                                        random_markers = TRUE,
                                        mutation_rate = 0.0,
                                        substitution_matrix =
                                          matrix(1 / 4, 4, 4)) {

  input_data_population_1 <- verify_genomeadmixr_data(input_data_population_1,
                                                      verbose = verbose)
  input_data_population_2 <- verify_genomeadmixr_data(input_data_population_2,
                                                      verbose = verbose)
  RcppParallel::setThreadOptions(num_threads)
  if (!is.na(critical_fst)) {
    stop_at_critical_fst <- TRUE
  }

  if (stop_at_critical_fst) {
    if (verbose) message("stopping at FST of: ", critical_fst)
    return(simulate_sequence_until(input_data_population_1 =
                                     input_data_population_1,
                                   input_data_population_2 =
                                     input_data_population_2,
                                   pop_size = pop_size,
                                   total_runtime = total_runtime,
                                   morgan = morgan,
                                   select_matrix = select_matrix,
                                   markers = markers,
                                   verbose = verbose,
                                   multiplicative_selection =
                                     multiplicative_selection,
                                   migration_rate = migration_rate,
                                   critical_fst = critical_fst,
                                   generations_between_update =
                                     generations_between_update,
                                   sampled_individuals = sampled_individuals,
                                   number_of_markers = number_of_markers,
                                   random_markers = random_markers,
                                   mutation_rate = mutation_rate,
                                   substitution_matrix = substitution_matrix,
                                   num_threads = num_threads,
                                   recombination_rate = recombination_rate))
  }

  if (!inherits(input_data_population_1, "genomeadmixr_data") ||
      !inherits(input_data_population_2, "genomeadmixr_data")) {
    stop("input_data should be of class genomeadmixr_data\n
          you can create such data with the functions\n
          create_input_data or vcfR_to_genomeadmixr_data")
  }

  select_matrix <- check_select_matrix(select_matrix,
                                       markers,
                                       use_data = TRUE)

  if (length(pop_size) == 1) {
    if (is.na(pop_size) ) {
      stop("need to provide population size")
    }
    pop_size <- c(pop_size, pop_size)
  }


  if (length(markers) == 1) {
    if (is.na(markers))  {
      markers <- c(-1, -1)
      track_frequency <- FALSE
    } else {
      markers <- c(markers)
      track_frequency <- TRUE
    }
  } else {
    track_frequency <- TRUE
  }

  if (track_frequency) {
    markers <- check_markers(markers, input_data_population_1$markers)
  }

  if (verbose) {
    message("markers: ", length(markers), "\n")
  }

  if (mutation_rate > 0) {
    verify_substitution_matrix(substitution_matrix)
  }

  recombination_map <- c(-1, -1)
  if (!is.na(recombination_rate)) {
    recombination_map <-
      create_recombination_map(input_data_population_1$markers,
                               recombination_rate)
  }

  selected_pop <- simulate_migration_emp_cpp(input_data_population_1$genomes,
                                             input_data_population_2$genomes,
                                             input_data_population_1$markers,
                                             select_matrix,
                                             pop_size,
                                             total_runtime,
                                             morgan,
                                             verbose,
                                             track_frequency,
                                             markers,
                                             multiplicative_selection,
                                             migration_rate,
                                             mutation_rate,
                                             substitution_matrix,
                                             num_threads,
                                             recombination_map)

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
                                         track_frequency)

  class(output) <- "genomadmixr_simulation"
  return(output)
}
