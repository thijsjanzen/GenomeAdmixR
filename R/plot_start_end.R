#' Plot both the starting frequencies and the final frequencies in one plot
#' @description This function plots the distribution of both the starting and
#' the final frequencies in one plot
#' @param results An object which is the result of \code{simulate_admixture},
#' being a list with four properties:
#' \code{population}, \code{frequencies}, \code{initial_frequencies} and
#' \code{final frequencies}
#' @param picked_ancestor Default is "ALL", where different colors indicate
#' different ancestors. Alternatively, for clarity, the user can specify a
#' specific ancestral allele, and only that allele is plotted
#' @param picked_population If multiple populations were simulated (in the case
#' of \code{simulate_admixture_migration}), which population should be plotted?
#' Default is population_1.
#' @return a ggplot object
#' @examples
#' markers <- seq(from = 0.2, to = 0.3, length.out = 100)
#'
#' pop <- simulate_admixture(pop_size = 1000,
#'                           number_of_founders = 3,
#'                           total_runtime = 11,
#'                           morgan = 1,
#'                           markers = markers)
#' require(ggplot2)
#' plot_start_end(pop,
#'                picked_ancestor = "ALL")
#' plot_start_end(pop,
#'                picked_ancestor = 1)
#' @export
plot_start_end <- function(results,
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

  a1_m <- dplyr::mutate(a1, timepoint = "start")
  a2_m <- dplyr::mutate(a2, timepoint = "end")

  to_plot_m <- rbind(a1_m, a2_m)

  if (picked_ancestor[[1]] == "ALL") {
    to_plot <- to_plot_m

    p1 <- ggplot2::ggplot(to_plot,
                        ggplot2::aes(x = .data[["location"]],
                                     y = .data[["frequency"]],
                                     colour = as.factor(.data[["ancestor"]]),
                                     group = interaction(.data[["ancestor"]],
                                                         .data[["timepoint"]]))) +
      ggplot2::geom_step(ggplot2::aes(lty = .data[["timepoint"]]))
  } else {

    to_plot <- subset(to_plot_m,
                      to_plot_m$ancestor %in% picked_ancestor)

    p1 <- ggplot2::ggplot(to_plot,
                          ggplot2::aes(x = .data[["location"]],
                                     y =  .data[["frequency"]],
                                     colour = as.factor(.data[["ancestor"]]),
                                     group = interaction(.data[["ancestor"]],
                                                         .data[["timepoint"]]))) +
      ggplot2::geom_step(ggplot2::aes(lty = .data[["timepoint"]]))
  }

  p1 <- p1 +
    ggplot2::xlab("Location (Morgan)") +
    ggplot2::ylab("Frequency") +
    ggplot2::labs(col = "Ancestor",
                  lty = "Time Point")

  return(p1)
}
