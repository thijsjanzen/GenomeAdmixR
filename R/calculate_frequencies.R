#' @keywords internal
calc_allele_frequencies <- function(indiv, alleles) {

  for (i in seq_along(indiv$chromosome1[, 1])) {
    left <- indiv$chromosome1[i, 1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome1[, 1])) {
      right <- indiv$chromosome1[i + 1, 1]
    }

    allele <- 1 + indiv$chromosome1[i, 2]
    alleles[allele] <- alleles[allele] + (right - left)
  }

  for (i in seq_along(indiv$chromosome2[, 1])) {
    left <- indiv$chromosome2[i, 1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome2[, 1])) {
      right <- indiv$chromosome2[i + 1, 1]
    }

    allele <- 1 + indiv$chromosome2[i, 2]
    alleles[allele] <- alleles[allele] + (right - left)
  }

  alleles <- alleles / sum(alleles)
  return(alleles)
}
