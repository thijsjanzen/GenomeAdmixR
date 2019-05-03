\name{calculate_LD}
\alias{calculate_LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate LD statistics
}
\description{
  This function calculates two matrices, once containing all pairwise LD values, one containing all pairwise r statistics
}
\usage{
calculate_LD(pop, sampled_individuals = 10,
                    number_of_markers = 100, random_markers = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop}{
  Focal population
}
  \item{sampled_individuals}{
  Number of individuals randomly sampled to calculate the LD matrices
}
  \item{number_of_markers}{
  Number of markers used to calculate the LD matrices
}
  \item{random_markers}{
  If TRUE, markers are randomly spaced along the chromosome, if FALSE, markers are equidistantly spaced along the chromosome.
}
}

\value{
 An object containing two items:
    \item{LD_matrix}{
       Pairwise LD statistics for all markers
    }
    \item{rsq_matrix}{
       Pairwise rsq statistics for all markers
    }
}

\author{
 Thijs Janzen
}

\examples{
  wildpop =  simulate_admixture(pop_size = 100,
                             number_of_founders = 10,
                             total_runtime = 100,
                             morgan = 1,
                             seed = 42)

  LD_results <- calculate_LD(wildpop,
                                    number_of_markers = 10,
                                    random_markers = TRUE)

  plot(LD_results$LD_matrix~LD_results$dist_matrix,
       pch = 16,
       xlab="Distance between markers",
       ylab = "Linkage Disequilibrium")
}
