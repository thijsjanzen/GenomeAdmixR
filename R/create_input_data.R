
#' @description Create data in a format that can be used by GenomeAdmixR
#' @param path path to input data
#' @param type type of data, options are 'ped', 'vcf', 'fasta'
#' @param chosen_chromosome GenomeAdmixR simulates only a single chromosome.
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