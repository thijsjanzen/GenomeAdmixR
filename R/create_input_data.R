#' read sequence data from file to be used in simulation
#' @description Create data in a format that can be used by GenomeAdmixR
#' @param file_names names of input files
#' @param type type of data, options are 'ped' and 'vcf'
#' @param chosen_chromosome GenomeAdmixR simulates only a single chromosome.
#' @param number_of_snps number of snps to be loaded from file, default
#' is to load all snps
#' @param random_snps if a subset of all snps has to be taken, should these
#' be sampled sequentially (e.g. the first 100 snps) or randomly (100 randomly
#' sampled snps) (examples are for 'number_of_snps' = 100).
#' @param verbose give verbose output
#' @return list with two properties: \code{genomes} a matrix with the
#' sequence translated to numerics, such that [actg] corresponds to [1234], and
#' missing data is represented with "-". Rows in the matrix correspond to
#' chromosomes, and columns represent bases. Two consecutive rows represent an
#' individual, such that rows 1-2 are individual, rows 3-4 are one individual
#' etc. \code{markers} corresponds to the locations of the markers (in bp) on
#' the chosen chromosome.
#' @export
read_input_data <- function(file_names,
                            type,
                            chosen_chromosome,
                            number_of_snps = NA,
                            random_snps = TRUE,
                            verbose = FALSE) {
  input_data <- c()
  if (type == "ped") {
    message("reading plink style data: ped/map pair")
    input_data <- read_ped(file_names[1], file_names[2], chosen_chromosome)
  }
  if (type == "vcf") {
    message("reading vcf data")
    input_data <- read_vcf(file_names[1], chosen_chromosome,
                           number_of_snps,
                           random_snps,
                           verbose)
  }
  return(input_data)
}

sample_output_matrix <- function(input_data_list,
                                 frequencies,
                                 pop_size) {

  num_markers <- length(input_data_list[[1]]$markers)

  output_matrix <- matrix(NA,
                          nrow = pop_size * 2,
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
  return(output_matrix)
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

  if (length(frequencies) == 1) {
    if (is.na(frequencies)) {
      frequencies <- rep(1 / length(input_data_list),
                         length(input_data_list))
    }
  }

  testit::assert(is.list(input_data_list))
  for (i in seq_along(input_data_list)) {
    testit::assert(inherits(input_data_list[[i]], "genomeadmixr_data"))
  }


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

  output <- list()
  output$genomes <- sample_output_matrix(input_data_list,
                                         frequencies,
                                         pop_size)
  output$markers <- input_data_list[[1]]$markers
  class(output) <- "genomeadmixr_data"
  return(output)
}

#' @keywords internal
convert_dna_to_numeric <- function(dna_matrix) {

  dna_matrix <- as.matrix(dna_matrix)
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

  odd_entries <- which(!(dna_matrix %in% c(0, 1, 2, 3, 4)))

  if (length(odd_entries) > 0) {
    message(paste0("found ", length(odd_entries), " non-biallelic entries"))
    message(paste0("these were set to missing data"))
    dna_matrix[odd_entries] <- 0
  }

  return(dna_matrix)
}

#' function to convert ped/map data to genome_admixr_data
#' @param simulation_data result of simulate_admixture
#' @param markers vector of locations of markers (in Morgan). If no vector is
#' provided, the function searches for marker locations in the simulation_data.
#' @param verbose provide verbose output (default is FALSE)
#' @return genomeadmixr_data object ready for simulate_admixture_data
#' @export
simulation_data_to_genomeadmixr_data <- function(simulation_data, # nolint
                                                 markers = NA,
                                                 verbose = FALSE) {

  output <- list()
  if (sum(is.na(markers)) > 0) {
    output$markers <- sort(unique(simulation_data$frequencies$location))
    if (is.null(output$markers)) {
      stop("no markers found, either provide them as argument, or make sure
           the simulation used to generate the data included markers")
    }
  } else {
    output$markers <- sort(unique(markers))
  }

  if (length(markers) > 1) {
    if (max(output$markers) <= 1.0) {
      # we have to rescale the markers to bp
      mark <- output$markers[output$markers > 0]
      min_mark <- min(mark)
      rescale_val <- 1 / min_mark
      output$markers <- output$markers * rescale_val
    }
  }

  if (!inherits(simulation_data, "population")) {
    if (is.list(simulation_data)) {
      if (inherits(simulation_data, "individual")) {
         indiv <- simulation_data
         simulation_data <- list()
         simulation_data[[1]] <- indiv
         class(simulation_data) <- "population"
      } else {
        if (inherits(simulation_data$population, "population")) {
          simulation_data <- simulation_data$population
        }
      }
    }
  }


  pop_for_cpp <- population_to_vector(simulation_data)
  gen_mat <- simulation_data_to_genomeadmixr_data_cpp(pop_for_cpp,
                                                      output$markers)

  output$genomes <- gen_mat
  if (sum(is.na(output$genomes)) > 0) {
    stop("NA values remain in genome matrix")
  }

  class(output) <- "genomeadmixr_data"
  return(output)
}

#' function to convert plink style (ped/map) data to genome_admixr_data
#' @param ped_data result of read.table(ped_file, header = F)
#' @param map_data result of read.table(map_file, header = F)
#' @param chosen_chromosome chromosome of choice
#' @param verbose verbose output
#' @return genomeadmixr_data object ready for simulate_admixture_data
#' @export
plink_to_genomeadmixr_data <- function(ped_data,  # nolint
                                       map_data,
                                       chosen_chromosome,
                                       verbose = FALSE) {
  # Base R
  map_data_duplicated <- map_data[rep(seq_len(nrow(map_data)), each = 2), ]

  map_indices <- which(map_data[, 1] == chosen_chromosome)
  map_data <- map_data[map_indices, ]

  ped_indices <- which(map_data_duplicated[, 1] == chosen_chromosome)
  ped_data <- ped_data[, 6 + ped_indices]

  # we have to first convert all letters to numbers,
  # then undouble
  if (verbose) message("converting atcg/ATCG entries to 1/2/3/4")
  ped_data <- convert_dna_to_numeric(ped_data)

  if (verbose) message("converting to numeric")
  ped_data <- as.matrix(ped_data,
                        nrow = length(ped_data[, 1]),
                        ncol = length(ped_data[1, ]))

  to_numeric <- function(v) {
    return(as.numeric(v))
  }
  ped_data <- apply(ped_data, 2, to_numeric)


  if (verbose) message("splitting markers into one row per chromosome")
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

  if (verbose) message("done")
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

  return(plink_to_genomeadmixr_data(ped_data,
                                    map_data,
                                    chosen_chromosome))
}

convert_vcf_to_alleles <- function(v) {
  if (is.na(v[[1]])) {
    return(c(0, 0))
  }

  available_alleles <- c(v[[2]], v[[3]])  # using ref / alt
  alleles <- stringr::str_split(v[[1]], pattern = "/")
  output <- c(NA, NA)
  for (i in 1:2) {
    output[i] <- available_alleles[1 + as.numeric(alleles[[1]][i])]
  }
  if (output[1] != output[2]) {
    # we don't know the phasing, so we shuffle
    if (stats::runif(1, 0, 1) < 0.5) {
      output <- rev(output)
    }
  }
  return(output)
}

convert_map_data_to_numeric <- function(map_data) {
  map_data <- map_data[, c(4, 5)]
  map_data[map_data == "a"] <- 1
  map_data[map_data == "A"] <- 1
  map_data[map_data == "c"] <- 2
  map_data[map_data == "C"] <- 2
  map_data[map_data == "t"] <- 3
  map_data[map_data == "T"] <- 3
  map_data[map_data == "TRUE"] <- 3
  map_data[map_data == "g"] <- 4
  map_data[map_data == "G"] <- 4
  map_data[is.na(map_data)] <- 0
  map_data[map_data == "0"] <- 0
  map_data[map_data == "?"] <- 0
  map_data[map_data == "N"] <- 0

  odd_entries <- which(!(map_data[, 1] %in% c("0", "1", "2", "3", "4")))
  odd_entries2 <- which(!(map_data[, 2] %in% c("0", "1", "2", "3", "4")))
  odd_entries <- unique(c(odd_entries, odd_entries2))

  map_data[odd_entries, 1] <- -1
  if (sum(map_data[, 1] == "-1") > 0) {
    message(paste0("found ", sum(map_data[, 1] == "-1"), " INDELs"))
    message(paste0("these were removed from the dataset"))
  }

  return(map_data)
}


#' function to convert a vcfR object to genome_admixr_data
#' @param vcfr_object result of vcfR::read.vcfR
#' @param chosen_chromosome chromosome of choice
#' @param number_of_snps number of snps to be loaded from the vcf file, default
#' is to load all snps
#' @param random_snps if a subset of all snps has to be taken, should these
#' be sampled sequentially (e.g. the first 100 snps) or randomly (100 randomly
#' sampled snps) (examples are for 'number_of_snps' = 100).
#' @param verbose if true, print progress bar
#' @return genomeadmixr_data object ready for simulate_admixture_data
#' @export
vcfR_to_genomeadmixr_data <- function(vcfr_object, chosen_chromosome,  # nolint
                                      number_of_snps = NA,
                                      random_snps = TRUE,
                                      verbose = FALSE) {

  if (!inherits(vcfr_object, "vcfR")) {
    stop("input has to be of class vcfR")
  }

  # now need to extract relevant data
  indices <- which(vcfr_object@fix[, 1] == chosen_chromosome)

  map_data <- convert_map_data_to_numeric(vcfr_object@fix[indices, ])
  to_remove <- which(map_data[, 1] == "-1")

  genome_data <- vcfr_object@gt[indices, ]
  genome_data <- genome_data[, -1]
  genome_data <- genome_data[-to_remove, ]

  marker_data <- vcfr_object@fix[indices, 2]
  marker_data <- marker_data[-to_remove]

  map_data <- map_data[-to_remove, ]
  v1 <- as.numeric(map_data[, 1])
  v2 <- as.numeric(map_data[, 2])
  map_data <- cbind(v1, v2)  # this is an ugly solution... but it works?

  if (!is.na(number_of_snps)) {
    # subsample snps
    snp_indices <- seq_along(map_data[, 1])
    if (random_snps) {
      snp_indices <- sort(sample(snp_indices, size = number_of_snps))
    } else {
      snp_indices <- 1:number_of_snps
    }
    map_data <- map_data[snp_indices, ]
    genome_data <- genome_data[snp_indices, ]
    marker_data <- marker_data[snp_indices]
  }

  if (verbose) message("extracting genotypes")

  numeric_matrix_for_cpp <- convert_to_numeric_matrix(genome_data)
  genome_matrix <- vcf_to_matrix_cpp(numeric_matrix_for_cpp,
                                     map_data[, 1],
                                     map_data[, 2])

  if (verbose) message("done")
  output <- list()
  output$genomes <- genome_matrix
  output$markers <- as.numeric(marker_data)
  class(output) <- "genomeadmixr_data"
  return(output)
}


#' @keywords internal
read_vcf <- function(vcf_name, chosen_chromosome,
                     number_of_snps,
                     random_snps, verbose = FALSE) {
  if (verbose) message("reading vcf file")
  vcf_data <- vcfR::read.vcfR(vcf_name)

  return(vcfR_to_genomeadmixr_data(vcf_data, chosen_chromosome,
                                   number_of_snps,
                                   random_snps, verbose))

}

#' function to generate artificial genomeadmixr_data
#' @param number_of_individuals number of individuals
#' @param marker_locations location of markers, either in bp or Morgan
#' @param used_nucleotides subset or full set of [1/2/3/4] (reflecting a/c/t/g)
#' @param nucleotide_frequencies frequencies of the used nucleotides, if not
#' provided, equal frequencies are assumed.
#' @return genomeadmixr_data object ready for simulate_admixture_data
#' @export
create_artificial_genomeadmixr_data <- function(number_of_individuals, # nolint
                                                marker_locations,
                                                used_nucleotides = 1:4,
                                                nucleotide_frequencies = NA) {

  if (is.na(nucleotide_frequencies)) {
    nucleotide_frequencies <- rep(1, length(used_nucleotides))
    nucleotide_frequencies <-
      nucleotide_frequencies / sum(nucleotide_frequencies)
  }
  num_markers <- length(marker_locations)
  fake_data <- list()
  fake_data$markers <- marker_locations

  if (length(used_nucleotides) == 1) {
    fake_data$genomes <- matrix(data = used_nucleotides,
                                nrow = number_of_individuals * 2,
                                ncol = num_markers)
  } else {
    fake_data$genomes <- matrix(data = sample(x = used_nucleotides,
                                              size = number_of_individuals * 2 *
                                                num_markers,
                                              replace = T),
                                nrow = number_of_individuals * 2,
                                ncol = num_markers)
  }

  class(fake_data) <- "genomeadmixr_data"
  return(fake_data)
}

#' function to write simulation output as PLINK style data
#' @param input_pop input population, either of class "population" or of class
#' "genomeadmixr_data"
#' @param marker_locations location of markers, in bp
#' @param file_name_prefix prefix of the ped/map files.
#' @param chromosome chromosome indication for map file
#' @param recombination_rate recombination rate in cM / kb
#' @return No return value
#' @export
write_as_plink <- function(input_pop,
                           marker_locations,
                           file_name_prefix,
                           chromosome = 1,
                           recombination_rate = 1) {


  if (is.list(input_pop)) {
    if (inherits(input_pop$population, "population")) {
      input_pop <- input_pop$population
    }
  }

  if (inherits(input_pop, "population")) {
    pop_for_cpp <- population_to_vector(input_pop)
    input_pop$genomes <- simulation_data_to_plink_cpp(pop_for_cpp,
                                                      marker_locations)
  } else {
    stop("input pop has to be of class population")
  }


  family_id <- "SIM"
  indiv_id <- paste0("indiv_", seq_along(input_pop$genomes[, 1]))
  paternal_id <- 0
  maternal_id <- 0
  sex  <- 0
  phenotype <- -9

  ped_table <- cbind(family_id,
                     indiv_id,
                     paternal_id,
                     maternal_id,
                     sex,
                     phenotype,
                     input_pop$genomes)

  utils::write.table(ped_table, file = paste0(file_name_prefix, ".ped"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  message("ped info written to: ", paste0(file_name_prefix, ".ped"))

  identifier <- paste0("rs", seq_along(marker_locations))
  recom_pos <- create_recombination_map(marker_locations, recombination_rate)
  recom_pos <- cumsum(recom_pos)

  output_matrix <- cbind(chromosome, identifier, recom_pos, marker_locations)
  utils::write.table(output_matrix, file = paste0(file_name_prefix, ".map"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  message("map info written to: ", paste0(file_name_prefix, ".map"))
}
