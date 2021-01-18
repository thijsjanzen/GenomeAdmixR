#' make a joy plot of the distribution of allele frequencies within a region
#' @description This function plots the distribution of allele frequencies
#' within a region over time, making use of a 'joyplot'
#' @param frequencies  A tibble containing four columns: \code{time},
#' \code{location}, \code{ancestor}, \code{frequency}. Typically one of the
#' items returned by \code{create_population_selection} or
#' \code{select_population} when the user specifies \code{track_frequency}.
#' @param time_points  A sequence of time points for which the user wants to
#' create the joyplot
#' @param picked_ancestor Default is "ALL", where different colors indicate
#' different ancestors. Alternatively, for clarity, the user can specify a
#' specific ancestral allele, and only that allele is plotted
#' @param picked_population If multiple populations were simulated (in the case
#' of \code{simulate_admixture_migration}), which population should be plotted?
#' Default is population_1.
#' @return a ggplot object
#' @export
#' @examples
#' \donttest{
#' s <- 0.01
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
#' plot_joyplot_frequencies(frequencies = selected_pop$frequencies,
#'                          time_points = 0:11,
#'                          picked_ancestor = "ALL")
#'
#' # joyplot frequencies returns a ggplot object, so we can
#' # add extra elements:
#' plot_joyplot_frequencies(frequencies = selected_pop$frequencies,
#'                          time_points = 0:11,
#'                          picked_ancestor = "ALL") +
#'   ggplot2::xlab("Location") +
#'   ggplot2::ylab("Generations")
#' }
#' @importFrom rlang .data
plot_joyplot_frequencies <- function(frequencies,
                                     time_points,
                                     picked_ancestor = "ALL",
                                     picked_population = 1
) {
  if ("population" %in% colnames(frequencies)) {
    frequencies <- subset(frequencies,
                          frequencies$population == picked_population)
  }

  time_points <- floor(time_points)
  vz <- subset(frequencies,
               frequencies$time %in% time_points)
  vz$ancestor <- as.factor(vz$ancestor)

  if (picked_ancestor == "ALL") {
    p1 <- ggplot2::ggplot(vz, ggplot2::aes(x = .data[["location"]],
                                           y = as.factor(.data[["time"]]),
                                           height = .data[["frequency"]],
                                           fill = .data[["ancestor"]])
                          ) +
      ggridges::geom_ridgeline(scale = 1.3)
  } else {
    vy <- subset(vz, vz$ancestor == picked_ancestor)
    p1 <- ggplot2::ggplot(vy, ggplot2::aes(x = .data[["location"]],
                                           y = as.factor(.data[["time"]]),
                                           height = .data[["frequency"]])) +
      ggridges::geom_ridgeline(scale = 1.3)
  }
  p1 <- p1 +
    ggplot2::labs(fill = "Ancestor")  +
    ggplot2::ylab("Time") +
    ggplot2::xlab("Location (Morgan)")
  return(p1)
}
