\name{calculate_tajima_d}
\alias{calculate_tajima_d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculate Tajima's d
}
\description{
 Tajima's d is calculated, given a number of markers.
}
\usage{
calculate_tajima_d(pop,
                   markers = seq(1e-6,1-1e-6,length.out = 100),
                   number_of_sampled_individuals = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop}{
  Population object
}
  \item{markers}{
  Vector from 0 to 1 (excluding 0 and 1) indicating the locations of the markers used for the analysis
}

  \item{number_of_sampled_individuals}{
  Number of individuals to base Tajima's d upon. Individuals are randomly drawn from the population.
}

}

\value{ A list with the following entries: \code{D} Tajima's D, \code{pi} pi, the average pairwise differences across the number of selected markers. \code{S} the number of segregating sites and \code{theta_hat}, the expected value of pi, calculated from S, such that theta = S/a1.
}

\author{
Thijs Janzen
}


\examples{
\dontrun{
pop <- simulate_admixture(pop_size = 100,
                          number_of_founders = 2,
                          total_runtime = 100,
                          seed = 42)$population

calculate_tajima_d(pop,
                   markers = seq(1e-6,1-1e-6,100),
                   number_of_sampled_individuals = 10)

}
}