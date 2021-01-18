#' Create isofemale
#' @description Creates isofemale individuals, given a population
#' @param source_pop Source population from which isofemales are generated
#' @param n Number of isofemales to be generated
#' @param inbreeding_pop_size Population size of the population used to generate
#' homozygous individuals
#' @param run_time Maximum runtime used for inbreeding
#' @param morgan Size of the chromosome in Morgan (e.g. the number of crossovers
#' during meiosis)
#' @param seed Random number generator seed
#' @param progress_bar Displays a progress_bar if TRUE. Default value is FALSE
#' @details To create an isofemale, two individuals are randomly picked from
#' the source population. Using these two individuals, a new population is
#' seeded, of size \code{inbreeding_pop_size}. Then, this population is allowed
#' to inbreed until either \code{run_time} is reached, or until all individuals
#' are homozygous and genetically identical, whatever happens first.
#' @return A list of length \code{n}, where each entry is a fully homozygous
#' isofemale.
#' @examples \donttest{
#' wildpop =  simulate_admixture(pop_size = 100,
#'                               number_of_founders = 10,
#'                               total_runtime = 5,
#'                               morgan = 1)
#'
#' isofemale <- create_iso_female(source_pop = wildpop,
#'                                n = 1,
#'                                inbreeding_pop_size = 100,
#'                                run_time = 100,
#'                                morgan = 1)
#' }
#' @export
create_iso_female <- function(source_pop,
                             n = 1,
                             inbreeding_pop_size = 100,
                             run_time = 2000,
                             morgan = 1,
                             seed = 42,
                             progress_bar = FALSE) {

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

    inbred_population <- simulate_admixture(input_population = parents,
                                            pop_size = inbreeding_pop_size,
                                            total_runtime = run_time,
                                            morgan = morgan)
    output_females[[i]] <-
          inbred_population$population[[
              sample(seq_along(inbred_population$population), 1)
                                      ]]

    class(output_females[[i]]) <- "individual"
  }
  return(output_females)
}
