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