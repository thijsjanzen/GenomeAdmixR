#' collect the full distribution of junctions in the population
#' @description calculates the distribution of junctions across the population
#' @param pop object of the class 'population'
#' @return vector with two entries per individual, each indicating the number of
#' junctions in the respective chromosomes
#' @export
calculate_dist_junctions <- function(pop) {
  get_num_junctions <- function(indiv) {
    v1 <- length(indiv$chromosome1[, 1]) - 1
    v2 <- length(indiv$chromosome2[, 1]) - 1 #subract one for start
    return(c(v1, v2))
  }

  vx <- unlist(lapply(pop, get_num_junctions))

  return(vx)
}

#' plot the distribution of junctions
#' @description plots the distribution of junctions in the population using
#' base R
#' @param pop of the class 'population'
#' @return No return value
#' @export
plot_dist_junctions <- function(pop) {
  junct <- calculate_dist_junctions(pop)
  vx <- table(junct)
  graphics::barplot(vx)
}
