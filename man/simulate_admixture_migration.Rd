\name{simulate_admixture_migration}
\alias{simulate_admixture_migration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Individual based simulation of the breakdown of contiguous ancestry blocks in two populations linked by migration
}
\description{
  Individual based simulation of the breakdown of contiguous ancestry blocks, with or without selection. Simulations can be started from scratch, or from a predefined input population. Two populations are simulated, connected by migration
}
\usage{
  simulate_admixture_migration(input_population_1 = NA,
                     input_population_2 = NA,
                     pop_size = c(100, 100),
                     initial_frequencies = list(c(1.0, 0),
                                                c(0, 1.0)),
                     total_runtime = 100,
                     morgan = 1,
                     seed,
                     select_matrix = NA,
                     markers = NA,
                     progress_bar = TRUE,
                     track_junctions = FALSE,
                     multiplicative_selection = TRUE,
                     migration_rate = 0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_population_1}{
    Potential earlier simulated population used as starting point for the simulation. If not provided by the user, the simulation starts from scratch.
  }
  \item{input_population_2}{
    Potential earlier simulated population used as starting point for the simulation. If not provided by the user, the simulation starts from scratch.
  }
  \item{pop_size}{
    Number of individuals in both populations.
  }
  \item{initial_frequencies}{
    A list describing the initial frequency of each ancestor in each population. Each entry in the list contains a vector with the frequencies for all ancestor. The length of the vector indicates the number of unique ancestors. If a vector not summing to 1 is provided, the vector is normalized.
  }
  \item{total_runtime}{
    Number of generations allowed for inbreeding
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

  \item{progress_bar}{
    Displays a progress_bar if TRUE. Default value is TRUE
  }

  \item{markers}{
    A vector of locations of markers (relative locations in [0, 1]). If a vector is provided, ancestry at these marker positions is tracked for every generation.
  }

  \item{track_junctions}{
    Track the average number of junctions over time if TRUE
  }

  \item{multiplicative_selection}{
    Default: TRUE. If TRUE, fitness is calculated for multiple markers by multiplying fitness values for each marker. If FALSE, fitness is calculated by adding fitness values for each marker.
  }
  \item{migration_rate}{
    Rate of migration between the two populations. Migration is implemented such that with probability m (migration rate) one of the two parents of a new offspring is from the other population, with probability 1-m both parents are of the focal population.
  }
}

\value{
  A list with: \code{population_1}, \code{population_2} two population objects, and three tibbles with allele frequencies (only contain values of a vector was provided to the argument \code{markers}: \code{frequencies} , \code{initial_frequencies} and \code{final_frequencies}. Each tibble contains five columns, \code{time}, \code{location}, \code{ancestor}, \code{frequency} and \code{population}, which indicates the number of generations, the location along the chromosome of the marker, the ancestral allele at that location in that generation, the frequency of that allele and the population in which it was recorded (1 or 2).
}

\author{
  Thijs Janzen
}

\examples{
\dontrun{
  select_matrix <- matrix(NA, nrow=1, ncol=5)

  s <- 0.1
  select_matrix[1, ] <- c(0.5, 0.5, 0.5+0.5*s, 0.5+s, 0)

  markers <- seq(from = 0.2, to = 0.6, by = 0.001)

  vy <- simulate_admixture_migration(seed = 42,
                                     migration_rate = 0.01,
                                     initial_frequencies = list(c(1,0),
                                                                c(0,1)),
                                     select_matrix = select_matrix,
                                     total_runtime = 100,
                                     markers = markers)
                                }
}