#' read sequence data from file to be used in simulation
#' @description Create data in a format that can be used by GenomeAdmixR
#' @param file_names names of input files
#' @param type type of data, options are 'ped' and 'vcf'
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
    message("reading vcf data")
    input_data <- read_vcf(file_names[1], chosen_chromosome)
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



  testit::assert(length(frequencies) == length(input_data_list))
  for (i in seq_along(input_data_list))  {
    for (j in seq_along(input_data_list)) {
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

  num_markers <- length(input_data_list[[1]]$markers)

  output_matrix <- matrix(NA, nrow = pop_size * 2,
                          ncol = num_markers)

  for (indiv in 1:pop_size) {
    chosen_pop <- sample(seq_along(input_data_list), size = 1,
                         prob = frequencies)

    focal_pop <- input_data_list[[chosen_pop]]$genomes
    indivs <- seq(from = 1, to = length(focal_pop[, 1]), by = 2)
    sampled_indiv <- sample(indivs, size = 1)
    indiv_index <- 1 + (indiv - 1) * 2
    output_matrix[indiv_index, ] <- focal_pop[sampled_indiv, ]
    output_matrix[indiv_index + 1, ] <- focal_pop[sampled_indiv + 1, ]
  }
  output <- list()
  output$genomes <- output_matrix
  output$markers <- input_data_list[[1]]$markers
  class(output) <- "genomeadmixr_data"
  return(output)
}

#' @keywords internal
convert_dna_to_numeric <- function(dna_matrix) {
  dna_matrix[dna_matrix == "a"] <- 1
  dna_matrix[dna_matrix == "A"] <- 1
  dna_matrix[dna_matrix == "c"] <- 2
  dna_matrix[dna_matrix == "C"] <- 2
  dna_matrix[dna_matrix == "t"] <- 3
  dna_matrix[dna_matrix == "T"] <- 3
  dna_matrix[dna_matrix == "TRUE"] <- 3
  dna_matrix[dna_matrix == "g"] <- 4
  dna_matrix[dna_matrix == "G"] <- 4
  dna_matrix[is.na(dna_matrix)] <- 0
  dna_matrix[dna_matrix == "0"] <- 0
  dna_matrix[dna_matrix == "?"] <- 0
  dna_matrix[dna_matrix == "N"] <- 0
  return(dna_matrix)
}

#' function to convert ped/map data to genome_admixr_data
#' @param ped_data result of read.table(ped_file, header = F)
#' @param map_data result of read.table(map_file, header = F)
#' @param chosen_chromosome chromosome of choice
#' @return genomeadmixr_data object ready for simulate_admixture_data
#' @export
ped_map_table_to_genomeadmixr_data <- function(ped_data,
                                               map_data,
                                               chosen_chromosome) {
  # Base R
  map_data_duplicated <- map_data[rep(seq_len(nrow(map_data)), each = 2), ]

  map_indices <- which(map_data[, 1] == chosen_chromosome)
  map_data <- map_data[map_indices, ]

  ped_indices <- which(map_data_duplicated[, 1] == chosen_chromosome)
  ped_data <- ped_data[, 6 + ped_indices]

  # we have to first convert all letters to numbers,
  # then undouble
  message("converting atcg/ATCG entries to 1/2/3/4")
  ped_data <- convert_dna_to_numeric(ped_data)

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

  for (j in seq_along(ped_data[, 1])) {
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

  class(output) <- "genomeadmixr_data"
  return(output)
}



#' @keywords internal
read_ped <- function(ped_name, map_name, chosen_chromosome) {
  message(paste("reading ped file:", ped_name))
  ped_data <- utils::read.table(ped_name, header = F)
  message(paste("reading map file:", map_name))
  map_data <- utils::read.table(map_name, header = F)

  return(ped_map_table_to_genomeadmixr_data(ped_data,
                                            map_data,
                                            chosen_chromosome))
}

convert_vcf_to_alleles <- function(v) {
  if (is.na(v[[1]])) {
    return(c(0, 0))
  }

  available_alleles <- c(v[[2]], v[[3]])  # ref / alt
  alleles <- stringr::str_split(v[[1]], pattern = "/")
  output <- c(NA, NA)
  for (i in 1:2) {
    output[i] <- available_alleles[ 1 + as.numeric(alleles[[1]][i]) ]
  }
  if (output[1] != output[2]) {
    if (stats::runif(1,0,1) < 0.5) { # we don't know the phasing, so we shuffle
      output <- rev(output)
    }
  }
  return(output)
}


#' function to convert a vcfR object to genome_admixr_data
#' @param vcfr_object result of vcfR::read.vcfR
#' @param chosen_chromosome chromosome of choice
#' @return genomeadmixr_data object ready for simulate_admixture_data
#' @export
vcfR_to_genomeadmixr_data <- function(vcfr_object, chosen_chromosome) {
  # now need to extract relevant data
  indices <- which(vcfr_object@fix[,1] == chosen_chromosome)

  map_data <- vcfr_object@fix[indices, ]
  genome_data <- vcfr_object@gt[indices, ]
  genome_data <- genome_data[, -1]

  num_indiv <- ncol(genome_data)
  marker_data <- as.numeric(map_data[, 2])
  num_markers <- length(marker_data)

  message("extracting genotypes")
  genome_matrix <- matrix(NA, nrow = num_indiv * 2, ncol = num_markers)
  for (i in seq_along(genome_data[1, ])) {
    indiv_seq <- genome_data[, i]
    to_convert <- cbind(indiv_seq, map_data[, 4], map_data[, 5])
    sequences <- apply(to_convert, 1, convert_vcf_to_alleles)
    indiv_index <- 1 + (i - 1) * 2
    genome_matrix[indiv_index, ] <- sequences[1, ]
    genome_matrix[indiv_index + 1, ] <- sequences[2, ]
  }

  message("converting atcg/ATCG entries to 1/2/3/4")
  genome_matrix <- convert_dna_to_numeric(genome_matrix)
  genome_matrix <- as.matrix(genome_matrix,
                             nrow = length(genome_matrix[, 1]),
                             ncol = length(genome_matrix[1, ]))

  to_numeric <- function(v) {
    return(as.numeric(v))
  }
  genome_matrix <- apply(genome_matrix, 2, to_numeric)
  message("done")
  output <- list()
  output$genomes <- genome_matrix
  output$markers <- marker_data
  class(output) <- "genomeadmixr_data"
  return(output)
}


#' @keywords internal
read_vcf <- function(vcf_name, chosen_chromosome) {
  message("reading vcf file")
  vcf_data <- vcfR::read.vcfR(vcf_name)

  return(vcfR_to_genomeadmixr_data(vcf_data, chosen_chromosome))

}