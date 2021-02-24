#' Create isofemale
#' @description Creates isofemale individuals, given a population
#' @param input_data Source population from which isofemales are generated
#' @param n Number of isofemales to be generated
#' @param inbreeding_pop_size Population size of the population used to generate
#' homozygous individuals
#' @param run_time Maximum runtime used for inbreeding
#' @param morgan Size of the chromosome in Morgan (e.g. the number of crossovers
#' during meiosis)
#' @param recombination_rate rate in cM / Mbp, used to map recombination to the
#' markers. If the recombination_rate is not set, the value for Morgan is used,
#' assuming that the markers included span an entire chromosome.
#' @param num_threads number of threads. Default is 1. Set to -1 to use all
#' available threads
#' @param verbose Displays verbose output if TRUE. Default value is FALSE
#' @details To create an isofemale, two individuals are randomly picked from
#' the source population. Using these two individuals, a new population is
#' seeded, of size \code{inbreeding_pop_size}. Then, this population is allowed
#' to inbreed until either \code{run_time} is reached, or until all individuals
#' are homozygous and genetically identical, whatever happens first.
#' @return A list of length \code{n}, where each entry is a fully homozygous
#' isofemale.
#' @export
iso_female_sequence <- function(input_data = NA,
                                n = 1,
                                inbreeding_pop_size = 100,
                                run_time = 2000,
                                morgan = 1,
                                recombination_rate = NA,
                                num_threads = 1,
                                verbose = FALSE) {

  input_data <- verify_genomeadmixr_data(input_data)

  indivs <- 1:(length(input_data$genomes[, 1]) /  2)

  # first we select the individuals that will be the parents of the isofemales
  indices <- sample(indivs,
                    n * 2,
                    replace = FALSE)

  # for each isofemale
  # pick the female, then simulate until run_time
  # then pick a random individual (assuming the population is fixed now)

  output_females <- list()
  for (i in 1:n) {
    parents <- list()

    index_indiv_1 <- 1 + (indices[i] - 1) * 2
    index_indiv_2 <- 1 + (indices[i + n] - 1) * 2


    parents$markers <- input_data$markers
    parents$genomes <- rbind(input_data$genomes[index_indiv_1, ],
                             input_data$genomes[index_indiv_1 + 1, ],
                             input_data$genomes[index_indiv_2, ],
                             input_data$genomes[index_indiv_2 + 1, ])
    class(parents) <- "genomeadmixr_data"

    inbred_population <- simulate_admixture(module = sequence_module(molecular_data = parents,
                                                                     morgan = morgan,
                                                                     recombination_rate = recombination_rate),
                                            pop_size = inbreeding_pop_size,
                                            total_runtime = run_time,
                                            verbose = verbose)
    output_females[[i]] <-
      inbred_population$population[[
        sample(seq_along(inbred_population$population), 1)
      ]]

    class(output_females[[i]]) <- "individual"
  }
  return(output_females)
}
