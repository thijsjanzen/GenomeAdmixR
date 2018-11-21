\name{create_iso_female}
\alias{create_iso_female}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Create isofemale
}
\description{
  Creates isofemale individuals, given a population
}
\usage{
create_iso_female(source_pop, n = 1,
                  inbreeding_pop_size = 100, run_time = 2000,
                  morgan = 1, seed = 42, progress_bar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_pop}{
  Source population from which isofemales are generated
}
  \item{n}{
  Number of isofemales to be generated
}
  \item{inbreeding_pop_size}{
  Population size of the population used to generate homozygous individuals
}
  \item{run_time}{
  Maximum runtime used for inbreeding
}
  \item{morgan}{
  Size of the chromosome in Morgan (e.g. the number of crossovers during meiosis)
}
  \item{seed}{
  Random number generator seed
  }
  \item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}
}
\details{
  To create an isofemale, two individuals are randomly picked from the source population. Using these two individuals, a new population is seeded, of size \code{inbreeding_pop_size}. Then, this population is allowed to inbreed until either \code{run_time} is reached, or until all individuals are homozygous and genetically identical, whatever happens first.
}
\value{
  A list of length \code{n}, where each entry is a fully homozygous isofemale.
}
\author{
 Thijs Janzen
}

\examples{
\dontrun{
  wildpop =  simulate_admixture(pop_size = 100,
                                number_of_founders = 10,
                                total_runtime = 5,
                                morgan = 1,
                                seed = 42,
                                progress_bar = TRUE)

  isofemale <- create_iso_female(source_pop = wildpop,
                               n = 1,
                               inbreeding_pop_size = 100,
                               run_time = 100,
                               morgan = 1,
                               seed = 42,
                               progress_bar = TRUE)
 }
}