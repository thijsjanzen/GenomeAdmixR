# #' check input population
# #' @param pop population, object to be checked
# #' @return corrected input population, or error if not supplied
#' @rawNamespace useDynLib(GenomeAdmixR, .registration = TRUE)
#' @rawNamespace import(Rcpp)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
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
check_select_matrix <- function(select_matrix,
                                markers,
                                use_data = FALSE) {
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

  if (dim(select_matrix)[2] == 5 && use_data == TRUE) {
    select_matrix[, 5] <- convert_dna_to_numeric(select_matrix[, 5])

    # this is super ugly code, but at least it works.
    other_matrix <- matrix(NA, nrow = length(select_matrix[, 1]),
                           ncol = 5)
    for (i in seq_along(select_matrix[, 1])) {
      for (j in 1:5) {
        other_matrix[i, j] <- as.numeric(select_matrix[i, j])
      }
    }
    select_matrix <- other_matrix

    sites_under_selection <- select_matrix[, 1]

    if (!(sites_under_selection %in% markers)) {
      stop("location of sites under selection have to exist in original data")
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

  get_chroms <- function(indiv) {
    a <- as.vector(t(indiv$chromosome1))
    b <- as.vector(t(indiv$chromosome2))
    return(c(a, b))
  }

  if (is.vector(source_pop)) return(source_pop)
  chroms <- lapply(source_pop, get_chroms)
  pop_for_cpp <- unlist(chroms)

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
create_random_markers <- function(number_of_markers,
                                  min_pos = 0,
                                  max_pos = 1) {
  markers <- c()
  while (length(markers) < number_of_markers) {
    temp_markers <- stats::runif(number_of_markers - length(markers),
                                 min_pos, max_pos)
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

  if (max(abs(indiv$chromosome1[, 2])) > 10000) {
    warning("Memory error recorded in chromosome\n")
    return(FALSE)
  }

  if (indiv$chromosome2[1, 1] != 0) {
    warning("Chromosome doesn't start at 0\n")
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
      warning("pop is nog of class population")
      return(FALSE)
    } else {
      return(verify_population(pop$population))
    }
  }
  v <- unlist(lapply(pop, verify_individual))
  if (sum(v) != length(v)) {
    warning("not all individuals passed test")
    return(FALSE)
  }

  return(TRUE)
}

# find local ancestry
# #' @description returns local ancestry at a site
# #' @param chrom chromosome
# #' @param pos position where to check for local ancestry
# #' @return local ancestry
#' @keywords internal
findtype <- function(chrom, pos) {

  if (pos < chrom[1, 1]) {
    return(NA)
  }

  if (pos == chrom[1, 1]) {
    return(chrom[1, 2])
  }

  if (pos >= utils::tail(chrom[, 1], 1)) {
    return(utils::tail(chrom[, 2], 1)) # return the last ancestry
  }

  b <- which(chrom[, 1] > pos)
  chromtype <- chrom[b[1] - 1, 2]

  return(chromtype[[1]])
}

#' print an individual to the console
#' @description prints an object of class genomeadmixr_data to the console
#' @param x individual
#' @param ... other arguments
#' @export
print.genomeadmixr_data <- function(x, ...) {
  print("Data to use as input for GenomeAdmixR")
  v1 <- paste("Number of individuals:",
              length(x$genomes[, 1]) / 2)
  v2 <- paste("Number of markers:",
              length(x$genomes[1, ]))
  print(v1)
  print(v2)
}


#' @keywords internal
print_substitution_matrix <- function(substitution_matrix) {
  for (i in 1:4) {
    print_str <- ""
    for (j in 1:4) {
      print_str <- paste(print_str, substitution_matrix[i, j])
    }
    message(print_str)
  }
}

#' @keywords internal
verify_substitution_matrix <- function(substitution_matrix) {
  if (length(substitution_matrix == 1)) {
    if (is.na(substitution_matrix[[1]])) {
      stop("No substitution matrix provided")
    }
  }

  if (dim(substitution_matrix)[[1]] != 4 ||
      dim(substitution_matrix)[[2]] != 4) {
    stop("\nCan not include mutations without proper substitution matrix\n",
         "substitution matrix should be a 4x4 matrix")
  }

  if (sum(substitution_matrix == 1) > 0) {
    warning("found rate matrix, rescaled all entries")
    # we have to rewrite as relative matrix
    for (i in 1:4) {
      row_entry <- substitution_matrix[i, ]
      row_entry <- row_entry * 0.25
      row_entry[i] <- 0
      row_entry[i] <- 1 - sum(row_entry)
      substitution_matrix[i, ] <- row_entry
    }
  }

  rs <- rowSums(substitution_matrix)
  if (sum(rs != 1) > 0) {
    warning("normalized rows to ensure they sum to 1")
    substitution_matrix[1:4, ] <- substitution_matrix[1:4, ] / rs[1:4]
  }

  message("using mutation with the following substitution matrix: ")
  print_substitution_matrix(substitution_matrix)

  return(substitution_matrix)
}

#' @keywords internal
check_markers <- function(markers, data_markers) {
  markers_in_data <- markers %in% data_markers
  which_markers_not_in_data <- which(markers_in_data == FALSE)

  if (length(which_markers_not_in_data) > 0) {
    warning(paste0("removing: ",
                   length(which_markers_not_in_data),
                   " markers as these do not exist in the original dataset"))
    markers <- markers[-which_markers_not_in_data]
  }

  if (length(markers) == 0) {
    stop("you have to provide markers that exist in the original dataset")
  }
  return(sort(markers))
}

#' @keywords internal
get_marker_range <- function(pop1, pop2) {
  get_min_pos <- function(indiv) {
    return(min(indiv$chromosome1[, 1],
               indiv$chromosome2[, 1]))
  }
  get_max_pos <- function(indiv) {
    return(max(indiv$chromosome1[, 1],
               indiv$chromosome2[, 1]))
  }

  min_positions1 <- unlist(lapply(pop1, get_min_pos))
  min_positions2 <- unlist(lapply(pop2, get_min_pos))

  max_positions1 <- unlist(lapply(pop1, get_max_pos))
  max_positions2 <- unlist(lapply(pop2, get_max_pos))

  min_pos <- min(min_positions1, min_positions2)
  max_pos <- max(max_positions1, max_positions2)
  return(c(min_pos, max_pos))
}

#' @keywords internal
check_for_bases <- function(pop) {

  unique_bases <- c()
  for (i in 1:10) {
    if (i < length(pop)) {
      b_c1 <- unique(pop[[i]]$chromosome1[, 2])
      b_c2 <- unique(pop[[i]]$chromosome2[, 2])
      unique_bases <- c(unique_bases, b_c1, b_c2)
    }
  }
  unique_bases <- sort(unique(unique_bases))
  using_sequencing_data <- FALSE
  # sum has to be at least 4, 0 represents missing data, but
  # might not be in the data
  if (sum(unique_bases %in% c(0, 1, 2, 3, 4)) >= 4 &&
      sum(unique_bases %in% c(0, 1, 2, 3, 4)) <= 5) {
    using_sequencing_data <- TRUE
  }
  return(using_sequencing_data)
}

#' @keywords internal
create_recombination_map <- function(markers,
                                     recombination_rate) {
  distances <- diff(markers)
  recomb_map <- c(0, distances)

  # recombination rate is in cM per Mbp
  rate_in_morgan <- recombination_rate / 100
  rate_per_bp <- rate_in_morgan / 1000000

  recomb_map <- recomb_map * rate_per_bp

  return(recomb_map)
}

#' @keywords internal
verify_genomeadmixr_data <- function(input_data, markers = NA) {
  if (!methods::is(input_data, "genomeadmixr_data")) {
    if (methods::is(input_data, "genomadmixr_simulation") ||
        methods::is(input_data, "individual")) {
      message("found simulation output, converting to genomeadmixr_data")
      message("this may take a while")
      input_data <-
        simulation_data_to_genomeadmixr_data(simulation_data =
                                               input_data,
                                             markers = markers)
      message("done converting, continuing as normal")
      return(input_data)
    } else {
      if (is.list(input_data)) {
        for (i in seq_along(input_data)) {
          input_data[[i]] <- verify_genomeadmixr_data(input_data[[i]], markers)
        }

        if (length(input_data) == 1) {
          input_data <- input_data[[1]]
        }

        return(input_data)
      }
    }
  }

  if (!methods::is(input_data, "genomeadmixr_data")) {
    input_data2 <- check_input_pop(input_data)
    if (methods::is(input_data2, "population")) {
      message("found simulation output, converting to genomeadmixr_data")
      message("this may take a while")
      input_data <-
        simulation_data_to_genomeadmixr_data(simulation_data = input_data,
                                             markers = markers)
      message("done converting, continuing as normal")
    } else {
      stop("input_data should be of class genomeadmixr_data
              you can create such data with the functions
              create_input_data or vcfR_to_genomeadmixr_data")
    }
  }
  return(input_data)
}

#' @keywords internal
convert_to_numeric_matrix <- function(genome_data) {

  genome_numeric <- genome_data
  genome_numeric[genome_numeric == "0/0"] <- 1
  genome_numeric[genome_numeric == "0/1"] <- 2
  genome_numeric[genome_numeric == "1/0"] <- 2
  genome_numeric[genome_numeric == "1/1"] <- 3
  genome_numeric[is.na(genome_numeric)]   <- 0

  genome_numeric <- t(genome_numeric)
  genome_numeric <- apply(genome_numeric, 2, as.numeric)

  return(genome_numeric)
}
