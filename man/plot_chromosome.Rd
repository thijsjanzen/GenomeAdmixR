% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_chromosome.R
\name{plot_chromosome}
\alias{plot_chromosome}
\title{plots a chromosome}
\usage{
plot_chromosome(chrom, xmin = 0, xmax = 1)
}
\arguments{
\item{chrom}{object of type chromosome, typically a table with two columns.
The first column indicates the start of an ancestry block (location in
Morgan), the second column indicates the ancestry type.}

\item{xmin}{minimum value of the range, default = 0.}

\item{xmax}{maximum value of the range, default = 1.}
}
\value{
No return value
}
\description{
This function plots a chromosome in the range [xmin, xmax].
Colors indicate different ancestry.
}
\examples{
wildpop =  simulate_admixture(
   module = ancestry_module(number_of_founders = 10, morgan = 1),
   pop_size = 1000,
   total_runtime = 10)

isofemale <- create_iso_female(
                 module = ancestry_module(input_population = wildpop,
                                          morgan = 1),
                 n = 1,
                 inbreeding_pop_size = 100,
                 run_time = 10)

plot_chromosome(chrom = isofemale[[1]]$chromosome1)
# and a detail of the chromosome:
plot_chromosome(chrom = isofemale[[1]]$chromosome1,
                xmin = 0.4,
                xmax = 0.6)
}
