# #' check input population
# #' @param pop population, object to be checked
# #' @return corrected input population, or error if not supplied
#' @rawNamespace useDynLib(GenomeAdmixR, .registration = TRUE)
#' @rawNamespace import(Rcpp)
#' @keywords internal
check_input_pop <- function(pop) {

  if (class(pop) == "individual") {
    pop <- list(pop)
    class(pop) <- "population"
  }

  if (!methods::is(pop, "population")) {
    if (is.list(pop)) {
      if (methods::is(pop$population, "population")) {
        pop <- pop$population
      } else{
        if (methods::is(pop$population_1, "population")) {
          warning("Warning, only population 1 is processed\n
              explicitly pass a population to remedy this\n")
          pop <- pop$population_1
        } else {
          types <- c()
          for (i in seq_along(pop)) {
            types[i] <- class(pop[[i]])
          }
          if (sum(types == "individual") == length(types)) {
            class(pop) <- "population"
          } else {
            stop("Input object is not of class 'population'")
          }
        }
      }
    } else {
      pop <- c(-1e6, -1e6)
    }
  }

  return(pop)
}

#' @keywords internal
check_initial_frequencies <- function(initial_frequencies) {
  if (!is.list(initial_frequencies)) {
    if (length(initial_frequencies) %% 2 == 1) {
      stop("wrong length of initial frequencies vector\n")
    }

    if (length(initial_frequencies) == 2) {
      warning("Assuming a single population with two ancestors\n")
      if (sum(initial_frequencies) != 1) {
        initial_frequencies <- initial_frequencies / sum(initial_frequencies)
        warning("initial frequencies were normalized to 1")
      }
    }

    message("found a vector instead of a list, converting by assuming")
    message("the second half is the second population\n")
    num_founders <- length(initial_frequencies) / 2
    output_freq <- list()
    a <- initial_frequencies[1:num_founders]
    output_freq[[1]] <- a
    b <- initial_frequencies[(num_founders + 1):length(initial_frequencies)]
    output_freq[[2]] <- b

    initial_frequencies <- output_freq
  }

  if (sum(initial_frequencies[[1]]) != 1) {
    initial_frequencies[[1]] <-
      initial_frequencies[[1]] / sum(initial_frequencies[[1]])
    warning("starting frequencies were normalized to 1\n")
  }
  if (sum(initial_frequencies[[2]]) != 1) {
    initial_frequencies[[2]] <-
      initial_frequencies[[2]] / sum(initial_frequencies[[2]])
    warning("starting frequencies were normalized to 1\n")
  }
  return(initial_frequencies)
}

#' @keywords internal
check_select_matrix <- function(select_matrix) {
  if (is.matrix(select_matrix)) {
    message(
      "Found a selection matrix, performing simulation including selection")
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
  return(select_matrix)
}

#' @keywords internal
generate_output_list_two_pop <- function(selected_pop,
                                         selected_popstruct_1,
                                         selected_popstruct_2,
                                         initial_freq_tibble,
                                         final_freq_tibble,
                                         track_frequency = FALSE,
                                         track_junctions = FALSE
) {
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
    colnames(selected_pop$frequencies) <- c("time",
                                            "location",
                                            "ancestor",
                                            "frequency",
                                            "population")
    frequencies_tibble <- tibble::as_tibble(selected_pop$frequencies,
                                            .name_repair = "unique")


    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if (track_frequency == TRUE && track_junctions == TRUE) {

    colnames(selected_pop$frequencies) <- c("time",
                                            "location",
                                            "ancestor",
                                            "frequency",
                                            "population")
    frequencies_tibble <- tibble::as_tibble(selected_pop$frequencies)

    output <- list("population_1" = selected_popstruct_1,
                   "population_2" = selected_popstruct_2,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }

  return(output)
}


#' @keywords internal
generate_output_list_one_pop <- function(selected_popstruct,
                                         selected_pop,
                                         initial_freq_tibble,
                                         final_freq_tibble,
                                         track_frequency = FALSE,
                                         track_junctions = FALSE) {
  output <- list()
  if (track_frequency == FALSE && track_junctions == FALSE) {
    output <- list("population" = selected_popstruct)
  }

  if (track_frequency == FALSE && track_junctions == TRUE) {
    output <- list("population" = selected_popstruct,
                   "junctions" = selected_pop$junctions)
  }

  if (track_frequency == TRUE && track_junctions == FALSE) {
    colnames(selected_pop$frequencies) <- c("time", "location",
                                            "ancestor", "frequency")
    frequencies_tibble <- tibble::as_tibble(selected_pop$frequencies)


    output <- list("population" = selected_popstruct,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if (track_frequency == TRUE && track_junctions == TRUE) {
    colnames(selected_pop$frequencies) <- c("time", "location",
                                            "ancestor", "frequency")

    frequencies_tibble <- tibble::as_tibble(selected_pop$frequencies)

    output <- list("population" = selected_popstruct,
                   "frequencies" = frequencies_tibble,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }
  return(output)
}



#' @keywords internal
population_to_vector <- function(source_pop) {
  if (is.vector(source_pop)) return(source_pop)
  pop_for_cpp <- c()
  for (i in seq_along(source_pop)) {
    x <- source_pop[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- source_pop[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }
  return(pop_for_cpp)
}

#' @keywords internal
increase_ancestor <- function(population, increment = 20) {
  increase_indiv <- function(indiv) {

    # -1 indicates the end, no increment there!
    positive <- which(indiv$chromosome1[, 2] > -1)
    indiv$chromosome1[positive, 2] <- indiv$chromosome1[positive, 2] + increment

    # -1 indicates the end, no increment there!
    positive <- which(indiv$chromosome2[, 2] > -1)
    indiv$chromosome2[positive, 2] <- indiv$chromosome2[positive, 2] + increment
    return(indiv)
  }

  if (!methods::is(population, "population")) {
    if (methods::is(population$population, "population")) {
      population <- population$population
    }
  }

  pop_2 <- lapply(population, increase_indiv)
  class(pop_2) <- "population"
  return(pop_2)
}

#' @keywords internal
create_random_markers <- function(number_of_markers) {
  markers <- c()
  while (length(markers) < number_of_markers) {
    temp_markers <- stats::runif(number_of_markers - length(markers), 0, 1)
    which_dupl <- which(duplicated(temp_markers))
    if (length(which_dupl)) {
      temp_markers <- temp_markers[-which_dupl]
    }
    markers <- c(markers, temp_markers)
  }
  markers <- sort(markers)
  return(markers)
}



#' @keywords internal
calc_allele_frequencies <- function(indiv, alleles) {

  for (i in seq_along(indiv$chromosome1[, 1])) {
    left <- indiv$chromosome1[i, 1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome1[, 1])) {
      right <- indiv$chromosome1[i + 1, 1]
    }

    allele <- 1 + indiv$chromosome1[i, 2]
    alleles[allele] <- alleles[allele] + (right - left)
  }

  for (i in seq_along(indiv$chromosome2[, 1])) {
    left <- indiv$chromosome2[i, 1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome2[, 1])) {
      right <- indiv$chromosome2[i + 1, 1]
    }

    allele <- 1 + indiv$chromosome2[i, 2]
    alleles[allele] <- alleles[allele] + (right - left)
  }

  alleles <- alleles / sum(alleles)
  return(alleles)
}

#' print an individual to the console
#' @description prints an object of class individual to the console
#' @param x individual
#' @param ... other arguments
#' @export
print.individual <- function(x, ...) {
  print("Individual with two Chromosomes")
  v1 <- paste("Chromosome 1:",
              length(x$chromosome1[, 1]) - 2,
              "junctions")
  v2 <- paste("Chromosome 2:",
              length(x$chromosome2[, 1]) - 2,
              "junctions")
  print(v1)
  print(v2)
}

#' print a population object
#' @description prints the contents of a population nicely
#' @param x input population
#' @param ... other arguments
#' @return No return value
#' @export
print.population <- function(x, ...) {
  v1 <- paste("Population with",
              length(x),
              "individuals")
  print(v1)
}

#' plot the genome of an individual
#' @description visualise ancestry blocks on both chromosomes
#' @param x object of type individual
#' @param cols colors for the different ancestors
#' @param ... other arguments
#' @return No return value
#' @export
plot.individual <- function(x, cols = NA,  ...) {
  alleles_chrom1 <- unique(x$chromosome1[, 2])
  alleles_chrom2 <- unique(x$chromosome2[, 2])
  num_colors <- 1 + max(alleles_chrom1, alleles_chrom2)
  if (num_colors > 20) num_colors <- 20
  color_palette <- grDevices::rainbow(num_colors, alpha = 1)

  if (!is.na(cols[[1]])) {
    color_palette <- cols
  }


  old_par <- graphics::par("mar", "mfrow", no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mfrow = c(2, 1))
  graphics::par(mar = c(2, 2, 2, 2))
  graphics::plot(NA,
                 xlim = c(0, 1),
                 ylim = c(0, 1),
                 xlab = "",
                 ylab = "",
                 xaxt = "n",
                 yaxt = "n",
                 bty  = "n")

  for (i in seq_along(x$chromosome1[, 1])) {
    xleft <- x$chromosome1[i, 1]
    xrght <- 1
    if (i < length(x$chromosome1[, 1])) {
      xrght <- x$chromosome1[i + 1, 1]
    }
    colour_index <- 1 + x$chromosome1[i, 2]
    colour_to_plot <- color_palette[colour_index]

    graphics::rect(xleft = xleft,
                   xright = xrght,
                   ybottom = 0,
                   ytop = 1,
                   col = colour_to_plot,
                   border = NA)
  }

  graphics::plot(NA,
                 xlim = c(0, 1),
                 ylim = c(0, 1),
                 xlab = "",
                 ylab = "",
                 xaxt = "n",
                 yaxt = "n",
                 bty  = "n")

  for (i in seq_along(x$chromosome2[, 1])) {
    xleft <- x$chromosome2[i, 1]
    xrght <- 1
    if (i < length(x$chromosome2[, 1])) {
      xrght <- x$chromosome2[i + 1, 1]
    }
    colour_index <- 1 + x$chromosome2[i, 2]
    colour_to_plot <- color_palette[colour_index]

    graphics::rect(xleft = xleft,
                   xright = xrght,
                   ybottom = 0,
                   ytop = 1,
                   col = colour_to_plot,
                   border = NA)
  }
}

#' @keywords internal
create_pop_class <- function(pop) {
  set_indiv_class <- function(indiv) {
    class(indiv) <- "individual"
    indiv
  }
  pop <- lapply(pop, set_indiv_class)
  class(pop) <- "population"
  return(pop)
}


# #' verify that an individual is correct
# #' @description function to verify correctness of an object of class
# #' 'individual'
# #' @param indiv object of class 'individual'
#' @keywords internal
verify_individual <- function(indiv) {

  if (!methods::is(indiv, "individual")) return(FALSE)

  if (indiv$chromosome1[1, 1] != 0) {
    warning("Chromosome doesn't start at 0\n")
    return(FALSE)
  }
  if (utils::tail(indiv$chromosome1, 1)[2] != -1) {
    warning("Chromosome doesn't end with -1\n")
    return(FALSE)
  }

  if (max(abs(indiv$chromosome1[, 2])) > 10000) {
    warning("Memory error recorded in chromosome\n")
    return(FALSE)
  }

  if (indiv$chromosome2[1, 1] != 0) {
    warning("Chromosome doesn't start at 0\n")
    return(FALSE)
  }

  if (utils::tail(indiv$chromosome2, 1)[2] != -1) {
    warning("Chromosome doesn't end with -1\n")
    return(FALSE)
  }

  if (max(abs(indiv$chromosome2[, 2])) > 10000) {
    warning("Memory error recorded in chromosome\n")
    return(FALSE)
  }

  return(TRUE)
}


# #' verify that a population is correct
# #' @description function to verify correctness of an object of class
# #' 'population'
# #' @param pop object of class 'population'
#' @keywords internal
verify_population <- function(pop) {

  if (!methods::is(pop, "population"))  {
    if (!methods::is(pop$population, "population")) {
      return(FALSE)
    } else {
      return(verify_population(pop$population))
    }
  }
  v <- unlist(lapply(pop, verify_individual))
  if (sum(v) != length(v)) return(FALSE)

  return(TRUE)
}

# find local ancestry
# #' @description returns local ancestry at a site
# #' @param chrom chromosome
# #' @param pos position where to check for local ancestry
# #' @return local ancestry
#' @keywords internal
findtype <- function(chrom, pos) {

  b <- which(chrom[, 1] > pos)
  chromtype <- chrom[b[1] - 1, 2]

  return(chromtype[[1]])
}
