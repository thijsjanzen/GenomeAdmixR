num_markers <- 100
num_indiv <- 100
chosen_markers <- 1:num_markers

fake_input_data1 <- list()
fake_input_data1$genomes <- matrix(data = 1,
                                   nrow = num_indiv,
                                   ncol = num_markers)


fake_input_data1$markers <- chosen_markers

fake_input_data2 <- list()
fake_input_data2$genomes <- matrix(data = 2,
                                   nrow = num_indiv,
                                   ncol = num_markers)
fake_input_data2$markers <- chosen_markers

class(fake_input_data1) <- "genomeadmixr_data"
class(fake_input_data2) <- "genomeadmixr_data"

select_matrix <- matrix(ncol = 5, nrow = 1)
s <- 1.0
select_matrix[1, ] <- c(50, 1.0, 1 + 0.5 * s, 1 + s, 1)

simul_two_pop <- simulate_admixture(
  module = sequence_module(molecular_data = list(fake_input_data1,
                                                 fake_input_data2),
                           markers = chosen_markers),
  migration = migration_settings(
    population_size = c(100, 100),
    migration_rate = 0.01),
  select_matrix = select_matrix,
  total_runtime = 100)