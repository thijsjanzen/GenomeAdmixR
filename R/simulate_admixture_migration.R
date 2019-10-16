#' Individual based simulation of the breakdown of contiguous ancestry blocks in
#' two populations linked by migration
#' @description IIndividual based simulation of the breakdown of contiguous
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
#' @param seed Seed of the pseudo-random number generator
#' @param select_matrix Selection matrix indicating the markers which are under
#' selection. If not provided by the user, the simulation proceeds neutrally.
#' If provided, each row in the matrix should contain five entries:
#' \code{location}{ location of the marker under selection (in Morgan) }
#' \code{fitness of wildtype (aa)} \code{fitness of heterozygote (aA)}
#' \code{fitness of homozygote mutant (AA)} \code{Ancestral type that
#' representes the mutant allele A}
#' @param progress_bar Displays a progress_bar if TRUE. Default value is TRUE
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
#' @return A list with: \code{population_1}, \code{population_2} two population
#' objects, and three tibbles with allele frequencies (only contain values of a
#' vector was provided to the argument \code{markers}: \code{frequencies},
#' \code{initial_frequencies} and \code{final_frequencies}. Each tibble contains
#' five columns, \code{time}, \code{location}, \code{ancestor}, \code{frequency}
#' and \code{population}, which indicates the number of generations, the
#' location along the chromosome of the marker, the ancestral allele at that
#' location in that generation, the frequency of that allele and the population
#' in which it was recorded (1 or 2).
#' @examples
#'  \dontrun{
#' select_matrix <- matrix(NA, nrow=1, ncol=5)
#'
#' s <- 0.1
#' select_matrix[1, ] <- c(0.5, 0.5, 0.5+0.5*s, 0.5+s, 0)
#'
#' markers <- seq(from = 0.2, to = 0.6, by = 0.001)
#'
#' vy <- simulate_admixture_migration(seed = 42,
#'                                   migration_rate = 0.01,
#'                                   initial_frequencies = list(c(1,0),
#'                                                              c(0,1)),
#'                                   select_matrix = select_matrix,
#'                                   total_runtime = 100,
#'                                   markers = markers)
#'}
#' @export
simulate_admixture_migration <- function(input_population_1 = NA,
                                         input_population_2 = NA,
                                         pop_size = c(100, 100),
                                         initial_frequencies = list(c(1.0, 0),
                                                                    c(0, 1.0)),
                                         total_runtime = 100,
                                         morgan = 1,
                                         seed = NULL,
                                         select_matrix = NA,
                                         markers = NA,
                                         progress_bar = TRUE,
                                         track_junctions = FALSE,
                                         multiplicative_selection = TRUE,
                                         migration_rate = 0.0) {

  cat("starting simulation incl migration\n")
  if (is.list(input_population_1)) {
    # if a list of individuals is given, the class is often wrong
    # let's check if that is the case
    if (!methods::is(input_population_1, "population")) {
      all_are_individuals <- vapply(input_population_1, class,
                                    FUN.VALUE = "character")
      if (sum(all_are_individuals == "individual") ==
         length(input_population_1)) {
        class(input_population_1) <- "population"
      }
    }

    if (methods::is(input_population_1$population, "population")) {
      input_population_1 <- input_population_1$population
    }

    if (methods::is(input_population_1, "population")) {
      input_population_1 <- population_to_vector(input_population_1)
    } else {
      input_population_1 <- c(-1e6, -1e6)
    }
  } else {
    input_population_1 <- c(-1e6, -1e6)
  }

  if (is.list(input_population_2)) {
    # if a list of individuals is given, the class is often wrong
    # let's check if that is the case
    if (!methods::is(input_population_2, "population")) {
      all_are_individuals <- vapply(input_population_2, class,
                                    FUN.VALUE = "character")
      if (sum(input_population_2 == "individual") ==
         length(input_population_2)) {
        class(input_population_2) <- "population"
      }
    }

    if (methods::is(input_population_2$population, "population")) {
      input_population_2 <- input_population_2$population
    }

    if (methods::is(input_population_2, "population")) {
      input_population_2 <- population_to_vector(input_population_2)
    } else {
      input_population_2 <- c(-1e6, -1e6)
    }
  } else {
    input_population_2 <- c(-1e6, -1e6)
  }

  if (sum(initial_frequencies[[1]]) != 1) {
    initial_frequencies[[1]] <-
      initial_frequencies[[1]] / sum(initial_frequencies[[1]])
    cat("starting frequencies were normalized to 1\n")
  }
  if (sum(initial_frequencies[[2]]) != 1) {
    initial_frequencies[[2]] <-
      initial_frequencies[[2]] / sum(initial_frequencies[[2]])
    cat("starting frequencies were normalized to 1\n")
  }

  if (is.matrix(select_matrix)) {
    cat("Found a selection matrix, performing simulation\n")
    cat("including selection\n")
    if (sum(is.na(select_matrix))) {
      stop("Can't start, there are NA values in the selection matrix!\n")
    }

    if (dim(select_matrix)[[2]] != 5) {
      stop("Incorrect dimensions of select_matrix,
           are you sure you provided all fitnesses?\n")
    }
    } else {
      if (is.na(select_matrix)) {
        select_matrix <- matrix(-1, nrow = 2, ncol = 2)
      }
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

  if (is.null(seed)) {
    seed <- round(as.numeric(Sys.time()))
  }

  selected_pop <- simulate_migration_cpp(input_population_1,
                                input_population_2,
                                select_matrix,
                                pop_size,
                                init_freq_matrix,
                                total_runtime,
                                morgan,
                                progress_bar,
                                track_frequency,
                                markers,
                                track_junctions,
                                multiplicative_selection,
                                migration_rate,
                                seed)

  selected_popstruct_1 <- create_pop_class(selected_pop$population_1)
  selected_popstruct_2 <- create_pop_class(selected_pop$population_2)

  initial_freq_tibble <- tibble::as.tibble(selected_pop$initial_frequencies)
  colnames(initial_freq_tibble) <- c("time",
                                     "location",
                                     "ancestor",
                                     "frequency",
                                     "population")

  final_freq_tibble <- tibble::as.tibble(selected_pop$final_frequencies)
  colnames(final_freq_tibble) <- c("time",
                                   "location",
                                   "ancestor",
                                   "frequency",
                                   "population")


  output <- list()
  if (track_frequency == FALSE && track_junctions == FALSE) {
    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2)
  }

  if (track_frequency == FALSE && track_junctions == TRUE) {
    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "junctions" = selected_pop$junctions)
  }

  if (track_frequency == TRUE && track_junctions == FALSE) {
    frequencies_tibble <- tibble::as.tibble(selected_pop$frequencies)
    colnames(frequencies_tibble) <- c("time",
                                      "location",
                                      "ancestor",
                                      "frequency",
                                        "population")

    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if (track_frequency == TRUE && track_junctions == TRUE) {
    frequencies_tibble <- tibble::as.tibble(selected_pop$frequencies)
    colnames(frequencies_tibble) <- c("time",
                                      "location",
                                      "ancestor",
                                      "frequency",
                                      "population")

    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }

  return(output)
}
