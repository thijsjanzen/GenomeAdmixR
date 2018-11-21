create_random_markers <- function(number_of_markers) {
  markers <- c()
  while (length(markers) < number_of_markers) {
    temp_markers <- runif(number_of_markers - length(markers), 0, 1)
    which_dupl <- which(duplicated(temp_markers))
    if (length(which_dupl)) {
      temp_markers <- temp_markers[-which_dupl]
    }
    markers <- c(markers, temp_markers)
  }
  markers <- sort(markers)
  return(markers)
}

calculate_dist_junctions <- function(pop) {
  get_num_junctions <- function(indiv) {
    v1 <- length(indiv$chromosome1[, 1]) - 2
    v2 <- length(indiv$chromosome2[, 1]) - 2 #subract one for start
    return(c(v1, v2))
  }

  vx <- unlist(lapply(pop, get_num_junctions))

  return(vx)
}

plot_dist_junctions <- function(pop) {
  junct <- calculate_dist_junctions(pop)
  vx <- table(junct)
  barplot(vx)
}

calculate_marker_frequency <- function(pop, location) {

  pop <- check_input_pop(pop)

  per_loc <- function(loc) {
    fun_chrom <- function(indiv) {
      return(c(findtype(indiv$chromosome1, loc),
               findtype(indiv$chromosome2, loc)))
    }
    types <- unlist(lapply(pop, fun_chrom))
    vv <- tibble::as.tibble(table(types))
    if(dim(vv)[1] > 0) {
      colnames(vv) <- c("ancestor", "frequency")
      vv$frequency <- vv$frequency / sum(vv$frequency)
      vv$location <- loc
      return(vv)
    } else {
      return(NA)
    }
  }

  all_types <- lapply(location, per_loc)
  output <- c()
  for(i in seq_along(all_types)) {
    output <- rbind(output, all_types[[i]])
  }
  output <- output[,c("location","ancestor","frequency")]

  return(output)
}

calculate_allele_frequencies <- function(source_pop,
                                         step_size,
                                         progress_bar = TRUE) {

  source_pop <- check_input_pop(source_pop)

  pop_for_cpp <- population_to_vector(source_pop)

  frequency_table <- calculate_allele_spectrum_cpp(pop_for_cpp,
                                                   step_size,
                                                   progress_bar)

  output <- tibble::as.tibble(frequency_table)
  colnames(output) <- c("location", "ancestor", "frequency")

  return(output)
}

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