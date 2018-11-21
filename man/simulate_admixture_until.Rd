\name{simulate_admixture_until}
\alias{simulate_admixture_until}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Simulate admixture until reaching threshold divergence
}
\description{
  Individual based simulation of the breakdown of contiguous ancestry blocks, with or without selection. Two populations are simulated, and upon userdefined timepoints, genetic divergence between the populations is checked. If the divergence exceeds a certain threshold, simulation is stopped.
}
\usage{
simulate_admixture_until(pop_size = 100,
                         number_of_founders = 2,
                         total_runtime = 100,
                         morgan = 1,
                         seed,
                         select_matrix = NA,
                         generations_between_update = 100,
                         critical_fst = 0.1,
                         sampled_individuals = 10,
                         number_of_markers = 100,
                         random_markers = TRUE,
                         overlap = 0,
                         multiplicative_selection = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop_size}{
  Number of individuals per population
}
  \item{number_of_founders}{
    Number of unique ancestors
}
  \item{total_runtime}{
  Maximum number of generations allowed to obtain level of genetic divergence
}
  \item{morgan}{
   Length of the chromosome in Morgan (e.g. the number of crossovers during meiosis)
}

\item{seed}{
  Seed of the pseudo-random number generator
}

 \item{select_matrix}{
  Selection matrix indicating the markers which are under selection. If not provided by the user, the simulation proceeds neutrally. If provided, each row in the matrix should contain five entries:
    \code{location}{
      location of the marker under selection (in Morgan)
    }
    \code{fitness of wildtype (aa)}
    \code{fitness of heterozygote (aA)}
    \code{fitness of homozygote mutant (AA)}
    \code{Ancestral type that representes the mutant allele A}
    }
  \item{generations_between_update}{
   Time interval used after which the FST is calculated again. Lower values increase accuracy, but also increase computation time
  }
  \item{critical_fst}{
   As long as the Fst between the two populations is lower than this fst, simulations continue
  }
  \item{sampled_individuals}{
    The number of individuals randomly sampled from the population to calculate the fst
  }
  \item{number_of_markers}{
    Number of markers used to calculate the Fst
  }
  \item{random_markers}{
    If TRUE, markers are randomly spaced along the chromosome, if FALSE, markers are equidistantly spaced along the chromosome.
  }
  \item{overlap}{
    Shared ancestry between the two populations. This is rounded to an integer, such that the number of shared ancestors is \code{round(number_of_founders * overlap)}
  }


\item{multiplicative_selection}{
  Default: TRUE. If TRUE, fitness is calculated for multiple markers by multiplying fitness values for each marker. If FALSE, fitness is calculated by adding fitness values for each marker.
}
}

\value{
A list with: \code{Population_1} a population object containing all individuals in population 1,\code{Population_2} a population object containing all individuals in population 2, \code{Number_of_generations} total number of generations required to obtain the cricital fst value, \code{FST} final FST value.
}

\author{
  Thijs Janzen
}
\examples{
  # should generate FST values around 0:
  example_pop <- simulate_admixture_until(pop_size = 100,
                                 number_of_founders = 2,
                                 total_runtime = 10,
                                 morgan = 1,
                                 seed = 42,
                                 generations_between_update = 1,
                                 critical_fst = 0.1,
                                 sampled_individuals = 10,
                                 number_of_markers = 100,
                                 random_markers = TRUE,
                                 overlap = 1,
                                 multiplicative_selection = TRUE)
}