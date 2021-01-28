#' Individual based simulation of the breakdown of contiguous ancestry blocks.
#' @description Individual based simulation of the breakdown of contiguous
#' ancestry blocks, with or without selection. Simulations can be started from
#' scratch, or from a predefined input population.
#' @param input_data Genomic data used as input, should be created by the
#' function \code{create_input_data} or by the function
#' \code{combine_input_data}
#' @param pop_size Vector containing the number of individuals in both
#' populations.
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
#' @param verbose Verbose output if TRUE. Default value is FALSE
#' @param markers A vector of locations of markers (relative locations in
#' [0, 1]). If a vector is provided, ancestry at these marker positions is
#' tracked for every generation.
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
#' @export
simulate_admixture_data <- function(input_data = NA,
                                    pop_size = NA,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = NULL,
                                    select_matrix = NA,
                                    markers = NA,
                                    verbose = FALSE,
                                    multiplicative_selection = TRUE) {

  if (class(input_data) != "genomeadmixr_data") {
    stop("input_data should be of class genomeadmixr_data\n
          you can create such data with the functions\n
          create_input_data or vcfR_to_genomeadmixr_data")
  }


  if (is.na(pop_size)) {
    if (length(input_data) == 2) {
      pop_size = length(input_data$genomes[, 1]) / 2
    } else {
      stop("pop_size is undefined, need an input population")
    }
  }

  select_matrix <- check_select_matrix(select_matrix)
  if (dim(select_matrix)[2] == 5) {
    select_matrix[, 5] <- convert_dna_to_numeric(select_matrix[, 5])

    # this is super ugly code, but at least it works.
    other_matrix <- matrix(NA, nrow = length(select_matrix[, 1]),
                           ncol = 5)
    for (i in 1:length(select_matrix[, 1])) {
      for (j in 1:5) {
        other_matrix[i, j] <- as.numeric(select_matrix[i, j])
      }
    }
    select_matrix <- other_matrix


    sites_under_selection <- select_matrix[, 1]
    normalized_markers <- (input_data$markers / max(input_data$markers)) * morgan
    if (!(sites_under_selection %in% normalized_markers)) {
      stop("location of sites under selection have to exist in original data")
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

  if (is.null(seed)) {
    seed <- round(as.numeric(Sys.time()))
  }


  selected_pop <- simulate_emp_cpp(input_data$genomes,
                                   input_data$markers,
                                   select_matrix,
                                   pop_size,
                                   total_runtime,
                                   morgan,
                                   verbose,
                                   track_frequency,
                                   markers,
                                   multiplicative_selection,
                                   seed)

  selected_popstruct <- create_pop_class(selected_pop$population)

  colnames(selected_pop$initial_frequencies) <- c("time",
                                                  "location",
                                                  "ancestor",
                                                  "frequency")
  colnames(selected_pop$final_frequencies) <-
    colnames(selected_pop$initial_frequencies)


  initial_freq_tibble <- tibble::as_tibble(selected_pop$initial_frequencies)
  final_freq_tibble   <- tibble::as_tibble(selected_pop$final_frequencies)


  output <- generate_output_list_one_pop(selected_popstruct,
                                         selected_pop,
                                         initial_freq_tibble,
                                         final_freq_tibble,
                                         track_frequency)

  output$frequencies <- convert_to_dna(output$frequencies)
  output$initial_frequency <- convert_to_dna(output$initial_frequency)
  output$final_frequency <- convert_to_dna(output$final_frequency)

  return(output)
}

#' @keywords internal
convert_to_dna <- function(frequency_tibble) {
  alleles <- frequency_tibble$ancestor
  alleles[alleles == 1] <- "a"
  alleles[alleles == 2] <- "c"
  alleles[alleles == 3] <- "t"
  alleles[alleles == 4] <- "g"
  alleles[alleles == 0] <- "-"

  frequency_tibble$ancestor <- alleles
  return(frequency_tibble)
}
