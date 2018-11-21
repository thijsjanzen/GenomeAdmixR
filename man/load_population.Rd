\name{load_population}
\alias{load_population}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Load a population from file
}
\description{
  Loads a population that has previously been written to file.
}
\usage{
load_population(file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{
      Name of the file that holds the population object
}
}
\details{
  This function is a wrapper for \code{readRDS}.
}
\value{
  A population object
}


\seealso{
\code{\link{save_population}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
wildpop <- simulate_admixture(pop_size = 10,
                               number_of_founders = 2,
                               total_runtime = 3,
                               morgan = 1,
                               seed = 123)

save_population(wildpop, file_name = "wildpop.pop")
wildpop2 <- load_population(file_name = "wildpop.pop")
all.equal(wildpop, wildpop2)
}