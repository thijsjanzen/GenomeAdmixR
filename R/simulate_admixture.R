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

  if(is.list(input_population)) {

    if(is(input_population$population, "population")) {
      input_population <- input_population$population
    }

    if(is(input_population, "population")) {
      input_population <- population_to_vector(input_population)
    } else {
      input_population <- c(-1e6, -1e6)
    }
  } else {
    input_population <- c(-1e6, -1e6)
  }

  if(sum(is.na(initial_frequencies))) {
    initial_frequencies <- rep(1.0 / number_of_founders,
                               times = number_of_founders)
  }

  if(sum(initial_frequencies) != 1) {
    initial_frequencies <- initial_frequencies / sum(initial_frequencies)
    cat("starting frequencies were normalized to 1\n")
  }

  if(is.matrix(select_matrix)) {
    if (sum(is.na(select_matrix))) {
      stop("Can't start, there are NA values in the selection matrix!\n")
    }

    if (dim(select_matrix)[[2]] != 5) {
      stop("Incorrect dimensions of select_matrix,
           are you sure you provided all fitnesses?\n")
    }
  } else {
    if(is.na(select_matrix)) {
      select_matrix <- matrix(-1, nrow=2,ncol=2)
    }
  }

  if(length(markers) == 1) {
    if(is.na(markers))  {
      markers <- c(-1, -1)
      track_frequency <- FALSE
    } else {
      track_frequency <- TRUE
    }
  } else {
    track_frequency <- TRUE
  }

  set.seed(seed)

  selected_pop <- simulate_cpp( input_population,
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
  colnames(initial_freq_tibble) <- c("location", "ancestor", "frequency")

  final_freq_tibble <- tibble::as.tibble(selected_pop$final_frequencies)
  colnames(final_freq_tibble) <- c("location", "ancestor", "frequency")


  output <- list()
  if(track_frequency == FALSE && track_junctions == FALSE) {
    output <- list("population" = selected_popstruct)
  }

  if(track_frequency == FALSE && track_junctions == TRUE) {
    output <- list("population" = selected_popstruct,
                   "junctions" = selected_pop$junctions)
  }

  if(track_frequency == TRUE && track_junctions == FALSE) {

    output <- list("population" = selected_popstruct,
                   "frequencies" = create_tibble_from_freq_table(
                                        selected_pop$frequencies,
                                        markers),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if(track_frequency == TRUE && track_junctions == TRUE) {

    output <- list("population" = selected_popstruct,
                   "frequencies" = create_tibble_from_freq_table(
                                        selected_pop$frequencies,
                                        markers),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }

  return(output)
}

create_admixed_individuals <- function(num_individuals,
                           population_size,
                           number_of_founders,
                           size_in_morgan) {

  pop <- create_pop_admixed_cpp(num_individuals,
                                number_of_founders,
                                population_size,
                                size_in_morgan)

  popstruct <- create_pop_class(pop$population)

  output <- list("population" = popstruct)
  return(output)
}



