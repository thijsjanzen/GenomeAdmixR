#' Function to manage settings associated with migration
#' @description creates a list with settings associated with migration.
#' @param migration_rate  Rate of migration between the two populations.
#' Migration is implemented such that with probability m (migration rate) one
#' of the two parents of a new offspring is from the other population, with
#' probability 1-m both parents are of the focal population.
#' @param stop_at_critical_fst option to stop at a critical FST value
#' , default is FALSE
#' @param critical_fst the critical fst value to stop, if
#' \code{stop_simulation_at_critical_fst} is TRUE
#' @param population_size vector of population sizes, one size for each
#' population
#' @param initial_frequencies A list describing the initial frequency of each
#' ancestor in each population. Each entry in the list contains a vector with
#' the frequencies for all ancestor. The length of the vector indicates the
#' number of unique ancestors. If a vector not summing to 1 is provided, the
#' vector is normalized.
#' @param generations_between_update The number of generations after which the
#' simulation has to check again whether the critical Fst value is exceeded
#' @param sampled_individuals Number of individuals to be sampled at random from
#' the population to estimate Fst
#' @param number_of_markers Number of markers to be used to estimate Fst
#' @param random_markers Are the markers to estimate Fst randomly distributed,
#' or regularly distributed? Default is TRUE.
#' @export
migration_settings <- function(migration_rate = NA,
                               stop_at_critical_fst = FALSE,
                               critical_fst = NA,
                               population_size = c(100, 100),
                               initial_frequencies = list(c(1.0, 0),
                                                          c(0, 1.0)),
                               generations_between_update = 10,
                               sampled_individuals = 10,
                               number_of_markers = 100,
                               random_markers = TRUE) {

  local_list = list(migration_rate = migration_rate,
                    stop_at_critical_fst = stop_at_critical_fst,
                    population_size = population_size,
                    initial_frequencies = initial_frequencies,
                    critical_fst = critical_fst,
                    generations_between_update = generations_between_update,
                    sampled_individuals = sampled_individuals,
                    number_of_markers = number_of_markers,
                    random_markers = random_markers,
                    type = "migration_settings")
  return(local_list)
}