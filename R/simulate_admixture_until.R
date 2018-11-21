simulate_admixture_until <- function(pop_size = 100,
                                     number_of_founders = 2,
                                     total_runtime = 100,
                                     morgan = 1,
                                     seed,
                                     select_matrix = NA,
                                     generations_between_update = 100,
                                     critical_fst = 0.1,
                                     sampled_individuals = 10,
                                     number_of_markers = 100,
                                     random_markers = TRUE,
                                     overlap = 0,
                                     multiplicative_selection = TRUE) {

  pop1 <- simulate_admixture(pop_size = pop_size,
                             number_of_founders = number_of_founders,
                             total_runtime = generations_between_update,
                             morgan = morgan,
                             seed = seed + 1,
                             select_matrix = select_matrix,
                             progress_bar = FALSE,
                             multiplicative_selection)$population

  pop2 <- simulate_admixture(pop_size = pop_size,
                             number_of_founders = number_of_founders,
                             total_runtime = generations_between_update,
                             morgan = morgan,
                             seed = seed + 2,
                             select_matrix = select_matrix,
                             progress_bar = FALSE,
                             multiplicative_selection)$population

  pop2 <- increase_ancestor(pop2,
                            increment = round(
                                  number_of_founders -
                                    number_of_founders * overlap))


  fst <- calculate_fst(pop1, pop2,
                       sampled_individuals = sampled_individuals,
                       number_of_markers = number_of_markers,
                       random_markers = TRUE)

  cnt <- 3
  cat("Number of Generations\tFST\n")
  cat(generations_between_update, "\t", fst, "\n")

  total_generations <- generations_between_update

  while(fst < critical_fst && total_generations < total_runtime) {
    pop1 <- simulate_admixture(pop1,
                               pop_size = pop_size,
                               total_runtime = generations_between_update,
                               morgan = morgan,
                               seed = seed + cnt,
                               select_matrix = select_matrix,
                               progress_bar = FALSE,
                               multiplicative_selection
                               )$population

    pop2 <- simulate_admixture(pop2,
                               pop_size = pop_size,
                               total_runtime = generations_between_update,
                               morgan = morgan,
                               seed = seed + cnt + 1,
                               select_matrix = select_matrix,
                               progress_bar = FALSE,
                               multiplicative_selection)$population
    cnt <- cnt + 2
    fst <- calculate_fst(pop1, pop2,
                         sampled_individuals = sampled_individuals,
                         number_of_markers = number_of_markers,
                         random_markers = TRUE)

    total_generations <- total_generations + generations_between_update
    cat(total_generations,"\t", fst ,"\n")
  }
  return(list("Population_1" = pop1,
              "Population_2" = pop2,
              "Number_of_generations" = total_generations,
              "FST" = fst))
}
