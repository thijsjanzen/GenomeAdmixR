#' read sequence data from file to be used in simulation
#' @description Create data in a format that can be used by GenomeAdmixR
#' @param file_names names of input files
#' @param type type of data, options are 'ped', 'vcf', 'fasta'
#' @param chosen_chromosome GenomeAdmixR simulates only a single chromosome.
#' @return list with two properties: \code{genomes} a matrix with the
#' sequence translated to numerics, such that [actg] corresponds to [1234], and
#' missing data is represented with "-". Rows in the matrix correspond to
#' chromosomes, and columns represent bases. Two consecutive rows represent an
#' individual, such that rows 1-2 are individual, rows 3-4 are one individual
#' etc. \code{markers} corresponds to the locations of the markers (in bp) on
#' the chosen chromosome.
#' @export
create_input_data <- function(file_names, type, chosen_chromosome) {
  input_data <- c()
  if (type == "ped") {
    message("reading plink style data: ped/map pair")
    input_data <- read_ped(file_names[1], file_names[2], chosen_chromosome)
  }
  if (type == "vcf") {
    mesage("reading vcf data")
    input_data <- read_vcf(file_names[1], chosen_chromosome)
  }
  if (type == "fasta") {
    message("reading fasta data")
    input_data <- read_fasta(file_names[1], chosen_chromosome)
  }
  return(input_data)
}

#' combine sequence data that was previously read from file into a population
#' @description Create data in a format that can be used by GenomeAdmixR,
#' entries are sampled randomly from each input data set, with replacement.
#' Probability of sampling from each input data set is driven by the input
#' frequencies, and total number of individuals sampled is driven by pop_size.
#' @param input_data_list list where each entry is the result of
#' \code{create_input_data}
#' @param frequencies frequency of each entry in the list in the starting
#' population
#' @param pop_size intended population size
#' @return the input data entries are combined to one single population that can
#' be used to seed \code{simulate_admixture_data}. Output is identical to
#' \code{create_input_data}
#' @export
combine_input_data <- function(input_data_list,
                               frequencies = NA,
                               pop_size) {

  output_matrix <- c()

  testit::assert(length(frequencies) == length(input_data_list))
  for (i in 1:length(input_data_list))  {
    for (j in 1:length(input_data_list)) {
      if (!all.equal(input_data_list[[i]]$markers,
                               input_data_list[[j]]$markers)) {
        stop("all input data sets need to use the same marker positions")
      }
    }
  }

  if (sum(frequencies) != 1) {
    message("frequencies normalized")
    frequencies <- frequencies / sum(frequencies)
  }


  for (indiv in 1:pop_size) {
    chosen_pop <- sample(seq_along(input_data_list), size = 1,
                         prob = frequencies)

    focal_pop <- input_data_list[[chosen_pop]]$genomes
    indivs <- seq(from = 1, to = length(focal_pop[, 1]), by = 2)
    sampled_indiv <- sample(indivs, size = 1)
    output_matrix <- rbind(output_matrix,
                           focal_pop[sampled_indiv, ],
                           focal_pop[sampled_indiv + 1, ])
  }
  output <- list()
  output$genomes <- output_matrix
  output$markers <- input_data_list[[1]]$markers
  return(output)
}


#' @keywords internal
read_ped <- function(ped_name, map_name, chosen_chromosome) {
  message(paste("reading ped file:", ped_name))
  ped_data <- read.table(ped_name, header = F)
  message(paste("reading map file:", map_name))
  map_data <- read.table(map_name, header = F)

  map_data_duplicated <- map_data[rep(seq_len(nrow(map_data)), each = 2), ]  # Base R


  map_indices <- which(map_data[, 1] == chosen_chromosome)
  map_data <- map_data[map_indices, ]

  ped_indices <- which(map_data_duplicated[, 1] == chosen_chromosome)
  ped_data <- ped_data[, 6 + ped_indices]

  # we have to first convert all letters to numbers,
  # then undouble
  message("converting atcg/ATCG entries to 1/2/3/4")
  ped_data[ped_data == "a"] <- 1
  ped_data[ped_data == "A"] <- 1
  ped_data[ped_data == "c"] <- 2
  ped_data[ped_data == "C"] <- 2
  ped_data[ped_data == "t"] <- 3
  ped_data[ped_data == "T"] <- 3
  ped_data[ped_data == "TRUE"] <- 3
  ped_data[ped_data == "g"] <- 4
  ped_data[ped_data == "G"] <- 4


  message("converting to numeric")
  ped_data <- as.matrix(ped_data,
                        nrow = length(ped_data[, 1]),
                        ncol = length(ped_data[1, ]))

  to_numeric <- function(v) {
    return(as.numeric(v))
  }
  ped_data <- apply(ped_data, 2, to_numeric)


  message("splitting markers into one row per chromosome")
  num_indiv <- length(ped_data[, 1])
  num_markers <- length(ped_data[1, ]) / 2

  output_matrix <- matrix(NA, ncol = num_markers,
                              nrow = num_indiv * 2)

  odd_indices <- seq(1, length(ped_data[1, ]), by = 2)
  even_indices <- seq(2, length(ped_data[1, ]), by = 2)

  for (j in 1:length(ped_data[, 1])) {
    chrom1 <- ped_data[j, odd_indices]
    chrom2 <- ped_data[j, even_indices]
    index1 <- 1 + (j - 1) * 2
    index2 <- index1 + 1
    output_matrix[index1, ] <- chrom1
    output_matrix[index2, ] <- chrom2
  }

  message("done")
  output <- list()
  output$genomes <- output_matrix
  output$markers <- map_data[, 4]
  return(output)
}