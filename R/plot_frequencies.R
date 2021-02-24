#' Plot the frequencies of all ancestors along the genome.
#' @description This function plots the frequency of all ancestors after
#' performing a simulation.
#' @param result An object which is the result of \code{select_population} or
#' \code{create_population_selection}, being a list with four properties:
#' \code{population}, \code{frequencies}, \code{initial_frequencies} and
#' \code{final frequencies}
#' @param locations A vector indicating the locations (in Morgan) where to
#' calculate the allele frequencies.
#' @param progress_bar Displays a progress_bar if TRUE. Default value is FALSE
#' @return a ggplot2 object
#' @examples
#' pop <- simulate_admixture(
#'              module = ancestry_module(number_of_founders = 4),
#'              pop_size = 1000,
#'              total_runtime = 11)
#' require(ggplot2)
#' plot_frequencies(result = pop)
#' @export
plot_frequencies <- function(result,
                             locations = seq(0, 1, length.out = 100),
                             progress_bar = FALSE) {

  to_plot <- calculate_allele_frequencies(result,
                                          locations,
                                          progress_bar)

  p1 <- ggplot2::ggplot(to_plot,
                        ggplot2::aes(x = .data[["location"]],
                                     y = .data[["frequency"]],
                                     colour = as.factor(.data[["ancestor"]]))) +
    ggplot2::geom_step()

  if (max(to_plot$location) < 10) {
    p1 <- p1 +
      ggplot2::xlab("Location (Morgan)")
  } else {
    p1 <- p1 +
      ggplot2::xlab("Location (bp)")
  }


  p1 <- p1 +
    ggplot2::ylab("Frequency") +
    ggplot2::labs(col = "Ancestor")

  return(p1)
}
