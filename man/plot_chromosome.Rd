\name{plot_chromosome}
\alias{plot_chromosome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 plot a chromosome
}
\description{
 This function plots a chromosome in the range [xmin, xmax]. Colors indicate different ancestry.
}
\usage{
plot_chromosome(chrom, xmin, xmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{
    object of type chromosome, typically a table with two columns. The first column indicates the start of an ancestry block (location in Morgan), the second column indicates the ancestry type.
}
  \item{xmin}{
    minimum value of the range, default = 0.
}
  \item{xmax}{
    maximum value of the range, default = 1.
}
}
\author{
  Thijs Janzen
}

\examples{
 wildpop =  simulate_admixture(pop_size = 100,
                             number_of_founders = 10,
                             total_runtime = 5,
                             morgan = 1,
                             seed = 42)

  isofemale <- create_iso_female(source_pop = wildpop,
                               n = 1,
                               inbreeding_pop_size = 100,
                               run_time = 10,
                               morgan = 1)

  plot_chromosome(isofemale[[1]]$chromosome1)
  plot_chromosome(isofemale[[1]]$chromosome1, xmin = 0.4, xmax = 0.6)
}

