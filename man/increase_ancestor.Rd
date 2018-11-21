\name{increase_ancestor}
\alias{increase_ancestor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 increase the ancestor index
}
\description{
 Increase the ancestor indexing to avoid overlap with another population
}
\usage{
increase_ancestor(population, increment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{population}{
    object of type 'chromosome population'
}
  \item{increment}{
    Number to add to the ancestor index, should ideally be larger than the number of ancestors in the other population. Default is 20, which should normally suffice.
}
}

\value{
 A population object
}

\author{
  Thijs Janzen
}

\examples{
wildpop <- simulate_admixture(pop_size = 100,
                             number_of_founders = 10,
                             total_runtime = 5,
                             morgan = 1,
                             seed = 42)

wildpop2 <- increase_ancestor(wildpop, 20)

}

