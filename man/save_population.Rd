\name{save_population}
\alias{save_population}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Save a population to file
}
\description{
  Saves a population to file for later use
}
\usage{
save_population(population, file_name, compression = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{population}{
    Object of class \code{population}
}
  \item{file_name}{
    Name of the file to save the population
}

  \item{compression}{
    By default, the population is compressed to reduce file size. See for more information \code{saveRDS}
}

}
\details{
    This function functions as a wrapper for the base function \code{saveRDS}.
}

\seealso{
\code{\link{load_population}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
wildpop <- simulate_admixture(pop_size = 10,
                               number_of_founders = 2,
                               total_runtime = 3,
                               morgan = 1,
                               seed = 123)

save_population(wildpop, file_name = "wildpop.pop")
}
