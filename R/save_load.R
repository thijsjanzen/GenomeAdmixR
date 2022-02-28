#' Save a population to file
#' @description Saves a population to file for later use
#' @param population Object of class \code{population}
#' @param file_name Name of the file to save the population
#' @param compression By default, the population is compressed to reduce file
#' size. See for more information \code{saveRDS}
#' @details This function functions as a wrapper for the base function
#' \code{saveRDS}.
#' @return No return value
#' @export
save_population <- function(population, file_name, compression = TRUE) {
  saveRDS(population, file = file_name, compress = compression)
}

#' Load a population from file
#' @description Loads a population that has previously been written to file.
#' @param file_name Name of the file to save the population
#' @return A population object
#' @seealso \code{\link{save_population}}
#' @details This function is a wrapper for \code{readRDS}.
#' @export
load_population <- function(file_name) {
  readRDS(file_name)
}
