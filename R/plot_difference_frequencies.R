#' Plot the change in frequency between the start and end of a simulation
#' @description This function plots the change in frequency of one or
#' multiple ancestors after performing a simulation.
#' @param results An object which is the result of \code{select_population} or
#' \code{create_population_selection}, being a list with four properties:
#' \code{population}, \code{frequencies}, \code{initial_frequencies} and
#' \code{final frequencies}
#' @param picked_ancestor Default is "ALL", where different colors indicate
#' different ancestors. Alternatively, for clarity, the user can specify a
#' specific ancestral allele, and only that allele is plotted
#' @param picked_population If multiple populations were simulated (in the case
#' of \code{simulate_admixture_migration}), which population should be plotted?
#' Default is population_1.
#' @return a ggplot2 object
#' @examples
#' s <- 0.1
#' select_matrix <- matrix(nrow = 1, ncol = 5)
#' select_matrix[1, ] <- c(0.25, 1.0, 1 + 0.5 * s, 1 + s, 0)
#'
#' markers <- seq(from = 0.2, to = 0.3, length.out = 100)
#'
#' selected_pop <- simulate_admixture(pop_size = 1000,
#'                                    number_of_founders = 10,
#'                                    total_runtime = 11,
#'                                    morgan = 1,
#'                                    select_matrix = select_matrix,
#'                                    markers = markers)
#' require(ggplot2)
#' plot_difference_frequencies(results = selected_pop,
#'                             picked_ancestor = "ALL")
#' @export
plot_difference_frequencies <- function(results,
                                        picked_ancestor = "ALL",
                                        picked_population = 1) {

  a1 <- results$initial_frequency
  a2 <- results$final_frequency

  if ("population" %in% colnames(a1)) {
    a1 <- subset(a1,
                 a1$population == picked_population)
    a2 <- subset(a2,
                 a2$population == picked_population)
  }

  a1 <- dplyr::select(a1, c("time", "location", "ancestor", "frequency"))
  a2 <- dplyr::select(a2, c("time", "location", "ancestor", "frequency"))

  colnames(a1) <- c("time",  "location", "ancestor", "frequency_before")
  colnames(a2) <- c("time",  "location", "ancestor", "frequency_after")

  ax <- dplyr::full_join(a1, a2, by = c("location", "ancestor"))

  ax_m <- dplyr::mutate(ax,
                        "diff_frequency" =
                          ax$frequency_after - ax$frequency_before)

  if (picked_ancestor[[1]] == "ALL") {
    to_plot <- ax_m

    p1 <- ggplot2::ggplot(to_plot,
                          ggplot2::aes(x = .data[["location"]],
                                       y = .data[["diff_frequency"]],
                                       colour = as.factor(.data[["ancestor"]]))) +
      ggplot2::geom_step()
  } else {

    to_plot <- subset(ax_m,
                      ax_m$ancestor %in% picked_ancestor)

    p1 <- ggplot2::ggplot(to_plot,
                          ggplot2::aes(x = .data[["location"]],
                                       y = .data[["diff_frequency"]],
                                       colour = as.factor(.data[["ancestor"]]))) +
      ggplot2::geom_step()
  }

  p1 <- p1 +
    ggplot2::xlab("Location (Morgan)") +
    ggplot2::ylab("Change in Frequency") +
    ggplot2::labs(col = "Ancestor")

  return(p1)
}
