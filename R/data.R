#' A subset of sequencing data from the Drosophila Genetics Reference Panel
#' @description This data set contains sequences from the 3R chromosome.
#' Included are 4603 SNPs with at least 0.05 minor allele frequency,  sequenced
#' across 410 isofemale lines. Sequences were downloaded from
#' <http://dgrp2.gnets.ncsu.edu/data.html>.
#' @references Mackay, T., Richards, S., Stone, E. et al. The Drosophila
#' melanogaster Genetic Reference Panel. Nature 482, 173–178 (2012).
#' <https://doi.org/10.1038/nature10811>
#' @format genomeadmixr_data object
#' @examples
#' data("dgrp2.3R.5k.data")
#' simulate_admixture(
#'        module = sequence_module(molecular_data = dgrp2.3R.5k.data),
#'        pop_size = 100,
#'        total_runtime = 10)
#' @usage data("dgrp2.3R.5k.data")
"dgrp2.3R.5k.data"
