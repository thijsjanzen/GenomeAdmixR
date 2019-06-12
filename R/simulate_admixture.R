#' Individual based simulation of the breakdown of contiguous ancestry blocks.
#' @description Individual based simulation of the breakdown of contiguous
#' ancestry blocks, with or without selection. Simulations can be started from
#' scratch, or from a predefined input population.
#' @param input_population Potential earlier simulated population used as
#' starting point for the simulation. If not provided by the user, the
#' simulation starts from scratch.
#' @param pop_size Vector containing the number of individuals in both
#' populations.
#' @param number_of_founders Number of unique ancestors
#' @param initial_frequencies A vector describing the initial frequency of each
#' ancestor. By default, equal frequencies are assumed. If a vector not summing
#' to 1 is provided, the vector is normalized.
#' @param total_runtime  Number of generations
#' @param morgan Length of the chromosome in Morgan (e.g. the number of
#' crossovers during meiosis)
#' @param seed Seed of the pseudo-random number generator
#' @param select_matrix Selection matrix indicating the markers which are under
#' selection. If not provided by the user, the simulation proceeds neutrally. If
#' provided, each row in the matrix should contain five entries:
#' \code{location}{ location of the marker under selection (in Morgan) }
#' \code{fitness of wildtype (aa)} \code{fitness of heterozygote (aA)}
#' \code{fitness of homozygote mutant (AA)} \code{Ancestral type that
#' represents the mutant allele A}
#' @param progress_bar Displays a progress_bar if TRUE. Default value is TRUE
#' @param markers A vector of locations of markers (relative locations in
#' [0, 1]). If a vector is provided, ancestry at these marker positions is
#' tracked for every generation.
#' @param track_junctions Track the average number of junctions over time if
#' TRUE
#' @param multiplicative_selection Default: TRUE. If TRUE, fitness is calculated
#' for multiple markers by multiplying fitness values for each marker. If FALSE,
#' fitness is calculated by adding fitness values for each marker.
#' @return A list with: \code{population} a population object, and three tibbles
#' with allele frequencies (only contain values of a vector was provided to the
#' argument \code{markers}: \code{frequencies} , \code{initial_frequencies} and
#' \code{final_frequencies}. Each tibble contains four columns, \code{time},
#' \code{location}, \code{ancestor} and \code{frequency}, which indicates the
#' number of generations, the location along the chromosome of the marker, the
#' ancestral allele at that location in that generation, and finally, the
#' frequency of that allele.
#' @examples
#' \dontrun{
#' wildpop <- simulate_admixture(pop_size = 10,
#'                              number_of_founders = 2,
#'                              total_runtime = 3,
#'                              morgan = 1,
#'                              seed = 123)
#'}
#' @export
simulate_admixture <- function(input_population = NA,
                               pop_size = 100,
                               number_of_founders = 2,
                               initial_frequencies = NA,
                               total_runtime = 100,
                               morgan = 1,
                               seed,
                               select_matrix = NA,
                               markers = NA,
                               progress_bar = TRUE,
                               track_junctions = FALSE,
                               multiplicative_selection = TRUE) {

  if (is.list(input_population)) {
    # if a list of individuals is given, the class is often wrong
    # let's check if that is the case
    if (!methods::is(input_population, "population")) {
      all_are_individuals <- vapply(input_population, class,
                                    FUN.VALUE = "character")
      if (sum(all_are_individuals == "individual") ==
            length(all_are_individuals)) {
        class(input_population) <- "population"
      }
    }

    if (methods::is(input_population$population, "population")) {
      input_population <- input_population$population
    }

    if (methods::is(input_population, "population")) {
      input_population <- population_to_vector(input_population)
    } else {
      input_population <- c(-1e6, -1e6)
    }
  } else {
    input_population <- c(-1e6, -1e6)
  }

  if (sum(is.na(initial_frequencies))) {
    initial_frequencies <- rep(1.0 / number_of_founders,
                               times = number_of_founders)
  }

  if (sum(initial_frequencies) != 1) {
    initial_frequencies <- initial_frequencies / sum(initial_frequencies)
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

  set.seed(seed)

  selected_pop <- simulate_cpp(input_population,
                               select_matrix,
                               pop_size,
                               number_of_founders,
                               initial_frequencies,
                               total_runtime,
                               morgan,
                               progress_bar,
                               track_frequency,
                               markers,
                               track_junctions,
                               multiplicative_selection)

  selected_popstruct <- create_pop_class(selected_pop$population)

  #initial_freq_tibble <- create_tibble_from_freq_mat(
  #                            selected_pop$initial_frequencies,
  #                            markers)

  #final_freq_tibble   <- create_tibble_from_freq_mat(
  #                            selected_pop$final_frequencies,
  #                            markers)

  initial_freq_tibble <- tibble::as.tibble(selected_pop$initial_frequencies)
  colnames(initial_freq_tibble) <- c("time",
                                     "location",
                                     "ancestor",
                                     "frequency")

  final_freq_tibble <- tibble::as.tibble(selected_pop$final_frequencies)
  colnames(final_freq_tibble) <- c("time", "location",
                                   "ancestor", "frequency")


  output <- list()
  if (track_frequency == FALSE && track_junctions == FALSE) {
    output <- list("population" = selected_popstruct)
  }

  if (track_frequency == FALSE && track_junctions == TRUE) {
    output <- list("population" = selected_popstruct,
                   "junctions" = selected_pop$junctions)
  }

  if (track_frequency == TRUE && track_junctions == FALSE) {
    frequencies_tibble <- tibble::as.tibble(selected_pop$frequencies)
    colnames(frequencies_tibble) <- c("time", "location",
                                      "ancestor", "frequency")

    output <- list("population" = selected_popstruct,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if (track_frequency == TRUE && track_junctions == TRUE) {
    frequencies_tibble <- tibble::as.tibble(selected_pop$frequencies)
    colnames(frequencies_tibble) <- c("time", "location",
                                      "ancestor", "frequency")

    output <- list("population" = selected_popstruct,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }

  return(output)
}
