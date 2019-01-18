\name{calculate_allele_frequencies}
\alias{calculate_allele_frequencies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculate allele frequencies
}
\description{
  Calculate for a number of regularly spaced markers the relative frequency of each ancestor in the population.
}
\usage{
calculate_allele_frequencies(source_pop, locations, progress_bar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_pop}{
  Population for which to estimate allele frequencies
}
  \item{locations}{
  A vector indicating the locations (in Morgan) where to calculate the allele frequencies.
}
\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}
}
\details{
  Markers are equidistantly spaced, with a distance of \code{step_size} in between them.
}
\value{
  A tibble containing the allele frequencies
}

\author{
  Thijs Janzen
}

\examples{
\dontrun{
  number_founders = 20
  wildpop =  simulate_admixture(pop_size = 1000,
                                number_of_founders = number_founders,
                                total_runtime = 10,
                                morgan = 1,
                                seed = 666)

 freq_output <- calculate_allele_frequencies(wildpop,
                                            step_size = 0.001,
                                            progress_bar = TRUE)

  require(ggplot2)
  ggplot(freq_output, aes(x=location, y = frequency, col = as.factor(ancestor))) +
    geom_line()
}
}