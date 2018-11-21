check_input_pop <- function(pop) {
  if(!is(pop, "population")) {
    if(is.list(pop)) {
      if(is(pop$population, "population")) {
        pop <- pop$population
      } else{
        stop("Input object is not of class 'population'")
      }
    } else {
      stop("Input object is not of class 'population'")
    }
  }
  return(pop)
}

population_to_vector <- function(source_pop) {
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

create_tibble_from_freq_table <- function(frequencies, select_matrix) {
  input_list <- list()
  for(i in 1:(dim(frequencies)[[3]])) {
    local_mat <- frequencies[ , , i]
    input_list[[i]] <- local_mat
  }

  to_apply <- function(local_mat, i) {
    time <- 0:(length(local_mat[[i]][,1])-1)
    marker_indicator <- rep(select_matrix[i], length(time))
    freq_tibble <- tibble::as.tibble( cbind(time,
                                            marker_indicator,
                                            local_mat[[i]]))
    colnames(freq_tibble) <- c("time",
                               "location",
                               0:(length(local_mat[[i]][1,]) - 1))

    freq_tibble <- tidyr::gather(freq_tibble,
                                 key = "ancestor",
                                 value = "frequency",
                                 -c(1,2))
    return(as.data.frame(freq_tibble))
  }

  interm_list <- lapply(seq_along(input_list), to_apply, local_mat = input_list)
  vy <- tibble::as.tibble(dplyr::bind_rows(interm_list))

  return(vy)
}

create_tibble_from_freq_mat <- function(frequencies, select_matrix) {
  found_markers <- c()
  for(i in 1:(dim(frequencies)[[1]])) {
    local_mat <- frequencies[i,]
    time <- 0
    marker_indicator <- select_matrix[i]
    freq_tibble <- c(time, marker_indicator, local_mat)

    found_markers <- rbind(found_markers, freq_tibble)
  }
  colnames(found_markers) <- c("time",
                               "location",
                               0:(length(frequencies[1,])-1))

  found_markers <- tibble::as.tibble(found_markers)
  found_markers <- tidyr::gather(found_markers,
                                 key = "ancestor",
                                 value = "frequency",
                                  -c(1,2))

  return(found_markers)
}

increase_ancestor <- function(population, increment = 20) {
  increase_indiv <- function(indiv) {

    # -1 indicates the end, no increment there!
    positive <- which(indiv$chromosome1[,2] > -1)
    indiv$chromosome1[positive, 2] <- indiv$chromosome1[positive, 2] + increment

    # -1 indicates the end, no increment there!
    positive <- which(indiv$chromosome2[,2] > -1)
    indiv$chromosome2[positive, 2] <- indiv$chromosome2[positive, 2] + increment
    return(indiv)
  }

  if(!is(population, "population")) {
    if(is(population$population, "population")) {
      population <- population$population
    }
  }

  pop_2 <- lapply(population, increase_indiv)
  class(pop_2) <- "population"
  return(pop_2)
}