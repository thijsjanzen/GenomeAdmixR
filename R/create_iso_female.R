#' function to simulate creation of an isofemale line
#' @description create_isofemale simulates the creation of an isofemale line
#' through extreme inbreeding.
#' @param module Source population from which isofemales are generated
#' @param n Number of isofemales to be generated
#' @param inbreeding_pop_size Population size of the population used to generate
#' homozygous individuals
#' @param run_time Maximum runtime used for inbreeding
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
create_iso_female <- function(module = ancestry_module(),
                             n = 1,
                             inbreeding_pop_size = 100,
                             run_time = 2000,
                             num_threads = 1,
                             verbose = FALSE) {
  if (module$type == "ancestry") {
    result <- iso_female_ancestry(source_pop = module$input_population,
                                  n = n,
                                  inbreeding_pop_size = inbreeding_pop_size,
                                  run_time = run_time,
                                  morgan = module$morgan,
                                  num_threads = num_threads,
                                  verbose = verbose)
    return(result)
  }
  if (module$type == "sequence") {
    result <- iso_female_sequence(input_data = module$input_data,
                                  n = n,
                                  inbreeding_pop_size = inbreeding_pop_size,
                                  run_time = run_time,
                                  morgan = module$morgan,
                                  recombination_rate = module$recombination_rate,
                                  num_threads = num_threads,
                                  verbose = verbose)
    return(result)
  }
  stop("wrong module used, please pick a module with
       ancestry_module() or sequence_module()")
}


#' Create isofemale
#' @description Creates isofemale individuals, given a population
#' @param source_pop Source population from which isofemales are generated
#' @param n Number of isofemales to be generated
#' @param inbreeding_pop_size Population size of the population used to generate
#' homozygous individuals
#' @param run_time Maximum runtime used for inbreeding
#' @param morgan Size of the chromosome in Morgan (e.g. the number of crossovers
#' during meiosis)
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
iso_female_ancestry <- function(source_pop = NA,
                              n = 1,
                              inbreeding_pop_size = 100,
                              run_time = 2000,
                              morgan = 1,
                              num_threads = 1,
                              verbose = FALSE) {

  source_pop <- check_input_pop(source_pop)

  # first we select the individuals that will be the parents of the isofemales
  indices <- sample(seq_along(source_pop), n * 2, replace = FALSE)

  #for each isofemale
  #pick the female, then simulate until run_time
  #then pick a random individual (assuming the population is fixed now)

  iso_females <- source_pop[indices]
  output_females <- list()
  for (i in 1:n) {
    parents <- list(iso_females[[i]], iso_females[[i + n]])
    class(parents) <- "population"

    inbred_population <- simulate_admixture(module = ancestry_module(
                                                      input_population = parents,
                                                      morgan = morgan
                                                    ),
                                            pop_size = inbreeding_pop_size,
                                            total_runtime = run_time,
                                            verbose = verbose,
                                            num_threads = num_threads)
    output_females[[i]] <-
          inbred_population$population[[
              sample(seq_along(inbred_population$population), 1)
                                      ]]

    class(output_females[[i]]) <- "individual"
  }
  return(output_females)
}
