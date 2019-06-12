#' Simulate admixture until
#' @description Individual based simulation of the breakdown of contiguous ancestry blocks, with or without selection. Two populations are simulated, and upon userdefined timepoints, genetic divergence between the populations is checked. If the divergence exceeds a certain threshold, simulation is stopped.
#' @param input_population_1 Potential earlier simulated population used as starting point for the simulation. If not provided by the user, the simulation starts from scratch.
#' @param input_population_2 Potential earlier simulated population used as starting point for the simulation. If not provided by the user, the simulation starts from scratch.
#' @param pop_size Vector containing the number of individuals in both populations.
#' @param initial_frequencies A list describing the initial frequency of each ancestor in each population. Each entry in the list contains a vector with the frequencies for all ancestor. The length of the vector indicates the number of unique ancestors. If a vector not summing to 1 is provided, the vector is normalized.
#' @param total_runtime  Number of generations
#' @param morgan Length of the chromosome in Morgan (e.g. the number of crossovers during meiosis)
#' @param seed Seed of the pseudo-random number generator
#' @param select_matrix Selection matrix indicating the markers which are under selection. If not provided by the user, the simulation proceeds neutrally. If provided, each row in the matrix should contain five entries: \code{location}{ location of the marker under selection (in Morgan) } \code{fitness of wildtype (aa)} \code{fitness of heterozygote (aA)} \code{fitness of homozygote mutant (AA)} \code{Ancestral type that representes the mutant allele A}
#' @param progress_bar Displays a progress_bar if TRUE. Default value is TRUE
#' @param markers A vector of locations of markers (relative locations in [0, 1]). If a vector is provided, ancestry at these marker positions is tracked for every generation.
#' @param track_junctions Track the average number of junctions over time if TRUE
#' @param multiplicative_selection Default: TRUE. If TRUE, fitness is calculated for multiple markers by multiplying fitness values for each marker. If FALSE, fitness is calculated by adding fitness values for each marker.
#' @param migration_rate  Rate of migration between the two populations. Migration is implemented such that with probability m (migration rate) one of the two parents of a new offspring is from the other population, with probability 1-m both parents are of the focal population.
#' @param critical_fst The threshold fst value after which the simulations need to stop.
#' @param generations_between_update The number of generations after which the simulation has to check again whether the critical Fst value is exceeded
#' @param sampled_individuals Number of individuals to be sampled at random from the population to estimate Fst
#' @param number_of_markers Number of markers to be used to estimate Fst
#' @param random_markers Are the markers to estimate Fst randomly distributed, or regularly distributed? Default is TRUE.
#' @return A list with: \code{Population_1} a population object containing all individuals in population 1,\code{Population_2} a population object containing all individuals in population 2, \code{Number_of_generations} total number of generations required to obtain the cricital fst value, \code{FST} final FST value.
#' @examples
#' \dontrun{
#'  should generate FST values around 0:
#' example_pop <- simulate_admixture_until(pop_size = c(100, 100),
#'                                        initial_frequencies = list(c(1.0, 1.0),
#'                                                                   c(1.0, 1.0)),
#'                                        total_runtime = 100,
#'                                        morgan = 1,
#'                                        seed = 42
#'                                        generations_between_update = 100,
#'                                        critical_fst = 0.1,
#'                                        sampled_individuals = 10,
#'                                        number_of_markers = 100,
#'                                        random_markers = TRUE)
#'}
#' @export
simulate_admixture_until <- function(input_population_1 = NA,
                                     input_population_2 = NA,
                                     pop_size = c(100, 100),
                                     initial_frequencies = list(c(1.0, 0),
                                                                c(0, 1.0)),
                                     total_runtime = 100,
                                     morgan = 1,
                                     seed,
                                     select_matrix = NA,
                                     markers = NA,
                                     progress_bar = TRUE,
                                     track_junctions = FALSE,
                                     multiplicative_selection = TRUE,
                                     migration_rate = 0.0,
                                     generations_between_update = 100,
                                     critical_fst = 0.1,
                                     sampled_individuals = 10,
                                     number_of_markers = 100,
                                     random_markers = TRUE) {

  pops <- simulate_admixture_migration(
    input_population_1 = input_population_1,
    input_population_2 = input_population_2,
    pop_size = pop_size,
    initial_frequencies = initial_frequencies,
    total_runtime = generations_between_update,
    seed = seed,
    morgan = morgan,
    select_matrix = select_matrix,
    markers = markers,
    progress_bar = progress_bar,
    track_junctions = track_junctions,
    multiplicative_selection = multiplicative_selection,
    migration_rate = migration_rate)

  fst <- calculate_fst(pops$population_1, pops$population_2,
                       sampled_individuals = sampled_individuals,
                       number_of_markers = number_of_markers,
                       random_markers = TRUE)

  cnt <- 3
  cat("Number of Generations\tFST\n")
  cat(generations_between_update, "\t", fst, "\n")

  total_generations <- generations_between_update

  while (fst < critical_fst && total_generations < total_runtime) {
    pops <- simulate_admixture_migration(
      input_population_1 = pops$population_1,
      input_population_2 = pops$population_2,
      pop_size = pop_size,
      initial_frequencies = initial_frequencies,
      total_runtime = generations_between_update,
      morgan = morgan,
      seed = seed + cnt,
      select_matrix = select_matrix,
      markers = markers,
      progress_bar = progress_bar,
      track_junctions = track_junctions,
      multiplicative_selection = multiplicative_selection,
      migration_rate = migration_rate)

    cnt <- cnt + 2
    fst <- calculate_fst(pops$population_1, pops$population_2,
                         sampled_individuals = sampled_individuals,
                         number_of_markers = number_of_markers,
                         random_markers = TRUE)

    total_generations <- total_generations + generations_between_update
    cat(total_generations, "\t", fst ,"\n")
  }
  return(list("Population_1" = pops$population_1,
              "Population_2" = pops$population_2,
              "Number_of_generations" = total_generations,
              "FST" = fst))
}
