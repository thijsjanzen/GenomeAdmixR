simulate_admixture_migration <- function(input_population_1 = NA,
                                         input_population_2 = NA,
                               pop_size = c(100, 100),
                               number_of_founders = 2,
                               initial_frequencies = list(NA, NA),
                               total_runtime = 100,
                               morgan = 1,
                               seed,
                               select_matrix = NA,
                               markers = NA,
                               progress_bar = TRUE,
                               track_junctions = FALSE,
                               multiplicative_selection = TRUE,
                               migration_rate = 0.0) {

  if(is.list(input_population_1)) {
    # if a list of individuals is given, the class is often wrong
    # let's check if that is the case
    if(!is(input_population_1, "population")) {
      all_are_individuals <- sapply(input_population_1, class)
      if(sum(input_population_1 == "individual") ==
         length(input_population_1)) {
        class(input_population_1) <- "population"
      }
    }

    if(is(input_population_1$population, "population")) {
      input_population_1 <- input_population_1$population
    }

    if(is(input_population_1, "population")) {
      input_population_1 <- population_to_vector(input_population)
    } else {
      input_population_1 <- c(-1e6, -1e6)
    }
  } else {
    input_population_1 <- c(-1e6, -1e6)
  }

  if(is.list(input_population_2)) {
    # if a list of individuals is given, the class is often wrong
    # let's check if that is the case
    if(!is(input_population_2, "population")) {
      all_are_individuals <- sapply(input_population_2, class)
      if(sum(input_population_2 == "individual") ==
         length(input_population_2)) {
        class(input_population_2) <- "population"
      }
    }

    if(is(input_population_2$population, "population")) {
      input_population_2 <- input_population_2$population
    }

    if(is(input_population_2, "population")) {
      input_population_2 <- population_to_vector(input_population)
    } else {
      input_population_2 <- c(-1e6, -1e6)
    }
  } else {
    input_population_2 <- c(-1e6, -1e6)
  }


  if(sum(is.na(initial_frequencies[[1]]))) {
    initial_frequencies[[1]] <- rep(1.0 / number_of_founders,
                               times = number_of_founders)
  }
  if(sum(is.na(initial_frequencies[[2]]))) {
    initial_frequencies[[2]] <- rep(1.0 / number_of_founders,
                                  times = number_of_founders)
  }

  if(sum(initial_frequencies[[1]]) != 1) {
    initial_frequencies[[1]] <- initial_frequencies[[1]] / sum(initial_frequencies)
    cat("starting frequencies were normalized to 1\n")
  }
  if(sum(initial_frequencies[[2]]) != 1) {
    initial_frequencies[[2]] <- initial_frequencies[[2]] / sum(initial_frequencies)
    cat("starting frequencies were normalized to 1\n")
  }


  if(is.matrix(select_matrix)) {
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

  init_freq <- c()
  for(i in 1:2) {
    local <- initial_frequencies[[i]]
    if(is.na(local)) {
      init_freq <- c(init_freq, c(0.5, 0.5))
    } else {
      init_freq <- c(init_freq, local)
    }
  }

  selected_pop <- simulate_migration_cpp( input_population_1,
                                input_population_2,
                                select_matrix,
                                pop_size,
                                number_of_founders,
                                init_freq,
                                total_runtime,
                                morgan,
                                progress_bar,
                                track_frequency,
                                markers,
                                track_junctions,
                                multiplicative_selection,
                                migration_rate)

  selected_popstruct_1 <- create_pop_class(selected_pop$population_1)
  selected_popstruct_2 <- create_pop_class(selected_pop$population_2)

  initial_freq_tibble <- tibble::as.tibble(selected_pop$initial_frequencies)
  colnames(initial_freq_tibble) <- c("time",
                                     "location",
                                     "ancestor",
                                     "frequency")

  final_freq_tibble <- tibble::as.tibble(selected_pop$final_frequencies)
  colnames(final_freq_tibble) <- c("time", "location", "ancestor", "frequency", "population")


  output <- list()
  if(track_frequency == FALSE && track_junctions == FALSE) {
    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2)
  }

  if(track_frequency == FALSE && track_junctions == TRUE) {
    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "junctions" = selected_pop$junctions)
  }

  if(track_frequency == TRUE && track_junctions == FALSE) {
    frequencies_tibble <- tibble::as.tibble(selected_pop$frequencies)
    colnames(frequencies_tibble) <- c("time","location","ancestor","frequency", "population")

    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if(track_frequency == TRUE && track_junctions == TRUE) {
    frequencies_tibble <- tibble::as.tibble(selected_pop$frequencies)
    colnames(frequencies_tibble) <- c("time","location","ancestor","frequency", "population")

    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }

  return(output)
}