\name{plot_frequencies}
\alias{plot_frequencies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plot the frequencies of all ancestors along the genome.
}
\description{
 This function plots the frequency of all ancestors after performing a simulation.
}
\usage{
plot_frequencies(results, locations, progress_bar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
    An object which is the result of \code{select_population} or \code{create_population_selection}, being a list with four properties: \code{population}, \code{frequencies}, \code{initial_frequencies} and \code{final frequencies}
}

  \item{locations}{
  A vector indicating the locations (in Morgan) where to calculate the allele frequencies.
}
\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is FALSE
}
}
\value{
  a ggplot2 object
}

\author{
  Thijs Janzen
}

\examples{
  s <- 0.1
  select_matrix <- matrix(nrow = 1, ncol = 5)
  select_matrix[1, ] <- c(0.5, 1.0, 1 + 0.5 * s, 1 + s, 0)

  selected_pop <- simulate_admixture(pop_size = 1000,
                                     number_of_founders = 10,
                                     total_runtime = 11,
                                     morgan = 1,
                                     select_matrix,
                                     seed = 1234)
  require(ggplot2)
  plot_frequencies(selected_pop, 0.01)
}

