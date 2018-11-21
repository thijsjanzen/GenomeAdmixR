\name{create_admixed_individuals}
\alias{create_admixed_individuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 create a population of admixed individuals
}
\description{
 Artificially create a population of admixed individuals
}
\usage{
create_admixed_individuals(num_individuals,
                           population_size,
                           number_of_founders,
                           size_in_morgan)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num_individuals}{
    the number of individuals to be created
}
  \item{population_size}{
    size of the population used to generate the individuals, this determines the degree of fixation (e.g. in smaller populations, individuals are more fixed for ancestry)
}
  \item{number_of_founders}{
    Number of unique ancestors
}
\item{size_in_morgan}{
   Length of the chromosome in Morgan (e.g. the number of crossovers during meiosis)
}
}

\value{
 A population object
}

\author{
  Thijs Janzen
}

\examples{
  admixed_pop <- create_admixed_individuals(num_individuals = 100,
                           population_size = 1e5,
                           number_of_founders = 2,
                           size_in_morgan = 1)
}