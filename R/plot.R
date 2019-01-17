joyplot_frequencies <- function(frequencies,
                                time_points,
                                picked_ancestor = "ALL"
                                )
{
  vz <- subset(frequencies,
               frequencies$time %in% time_points)
  vz$ancestor <- as.factor(vz$ancestor)

  if(picked_ancestor == "ALL") {
    p1 <- ggplot2::ggplot(vz, ggplot2::aes(x = vz$location,
                         y = as.factor(vz$time),
                         height = vz$frequency,
                         fill = vz$ancestor)) +
            ggridges::geom_ridgeline(scale = 1.3) +
            ggplot2::ylab("Time")
    return(p1)
  } else {
    vy <- subset(vz, vz$ancestor == picked_ancestor)
    p1 <- ggplot2::ggplot(vy, ggplot2::aes(x = vy$location,
                         y = as.factor(vy$time),
                         height = vy$frequency)) +
            ggridges::geom_ridgeline(scale = 1.3,
                                     fill = "lightblue") +
            ggplot2::ylab("Time")
    return(p1)
  }
}

plot_start_end <- function(results,
                           picked_ancestor = "ALL") {

  a1 <- results$initial_frequency
  a2 <- results$final_frequency

  a1_m <- dplyr::mutate(a1, timepoint = "start")
  a2_m <- dplyr::mutate(a2, timepoint = "end")

  to_plot_m <- rbind(a1_m, a2_m)

  if(picked_ancestor[[1]] == "ALL") {
    to_plot <- to_plot_m

    p1 <- ggplot2::ggplot(to_plot, ggplot2::aes(x = to_plot$location,
                              y = to_plot$frequency,
                              colour = to_plot$ancestor,
                              group = interaction(to_plot$ancestor,
                                                  to_plot$timepoint))) +
      ggplot2::geom_step(ggplot2::aes(lty = to_plot$timepoint))
  } else {

    to_plot <- dplyr::filter(to_plot_m,
                             to_plot_m$ancestor %in% picked_ancestor)

    p1 <- ggplot2::ggplot(to_plot, ggplot2::aes(x = to_plot$location,
                              y = to_plot$frequency,
                              colour = to_plot$ancestor,
                              group = interaction(to_plot$ancestor,
                                                  to_plot$timepoint))) +
      ggplot2::geom_step(ggplot2::aes(lty = to_plot$timepoint))
  }

  p1 <- p1 +
    ggplot2::xlab("Location (Morgan)") +
    ggplot2::ylab("Frequency") +
    ggplot2::labs(col = "Ancestor",
         lty = "Time Point")


  return(p1)
}

plot_difference_frequencies <- function(results,
                                        picked_ancestor = "ALL") {

  a1 <- results$initial_frequency
  a1 <- dplyr::select(a1, c("time","location","ancestor","frequency"))
  a2 <- results$final_frequency
  a2 <- dplyr::select(a2, c("time","location","ancestor","frequency"))

  colnames(a1) <- c("time"    ,  "location" , "ancestor" , "frequency_before")
  colnames(a2) <- c("time"    ,  "location" , "ancestor" , "frequency_after")

  ax <- dplyr::full_join(a1,a2, by = c("time", "location", "ancestor"))

  ax_m <- dplyr::mutate(ax,
                        "diff_frequency" =
                          ax$frequency_after - ax$frequency_before)

  if(picked_ancestor[[1]] == "ALL") {
    to_plot <- ax_m

    p1 <- ggplot2::ggplot(to_plot, ggplot2::aes(x = to_plot$location,
                                                y = to_plot$diff_frequency,
                                                colour = to_plot$ancestor)) +
      ggplot2::geom_step()
  } else {

    to_plot <- dplyr::filter(ax_m,
                             ax_m$ancestor %in% picked_ancestor)

    p1 <- ggplot2::ggplot(to_plot, ggplot2::aes(x = to_plot$location,
                                                y = to_plot$diff_frequency,
                                                colour = to_plot$ancestor)) +
      ggplot2::geom_step()
  }

  p1 <- p1 +
    ggplot2::xlab("Location (Morgan)") +
    ggplot2::ylab("Change in Frequency") +
    ggplot2::labs(col = "Ancestor")

  return(p1)
}

plot_frequencies <- function(result,
                             locations = seq(0, 1, length.out = 100),
                             progress_bar = FALSE) {

  to_plot <- calculate_allele_frequencies(result,
                                          locations,
                                          progress_bar)

  p1 <- ggplot2::ggplot(to_plot,
                          ggplot2::aes(x = to_plot$location,
                                       y = to_plot$frequency,
                                       colour = as.factor(to_plot$ancestor))) +
        geom_step()

  p1 <- p1 +
    ggplot2::xlab("Location (Morgan)") +
    ggplot2::ylab("Frequency") +
    ggplot2::labs(col = "Ancestor")

  return(p1)
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
