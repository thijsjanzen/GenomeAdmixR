#' plots a chromosome
#' @description This function plots a chromosome in the range [xmin, xmax].
#' Colors indicate different ancestry.
#' @param chrom  object of type chromosome, typically a table with two columns.
#' The first column indicates the start of an ancestry block (location in
#' Morgan), the second column indicates the ancestry type.
#' @param xmin minimum value of the range, default = 0.
#' @param xmax maximum value of the range, default = 1.
#' @return No return value
#' @examples
#' wildpop =  simulate_admixture(
#'    module = ancestry_module(number_of_founders = 10, morgan = 1),
#'    pop_size = 1000,
#'    total_runtime = 10)
#'
#' isofemale <- create_iso_female(
#'                  module = ancestry_module(input_population = wildpop,
#'                                           morgan = 1),
#'                  n = 1,
#'                  inbreeding_pop_size = 100,
#'                  run_time = 10)
#'
#' plot_chromosome(chrom = isofemale[[1]]$chromosome1)
#' # and a detail of the chromosome:
#' plot_chromosome(chrom = isofemale[[1]]$chromosome1,
#'                 xmin = 0.4,
#'                 xmax = 0.6)
#' @export
plot_chromosome <- function(chrom, xmin = 0, xmax = 1) {
  alleles <- unique(chrom[, 2])
  num_colors <- 1 + max(alleles)
  if (num_colors > 20) num_colors <- 20
  color_palette <- grDevices::rainbow(num_colors, alpha = 1)

  if (max(chrom[, 1]) > 1 && xmax == 1) {
    xmax <- max(chrom[, 1])
  }

  graphics::plot(NA,
       xlim = c(xmin, xmax),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")

  for (i in seq_along(chrom[, 1])) {
    xleft <- chrom[i, 1]
    xrght <- 1
    if (i < length(chrom[, 1])) {
      xrght <- chrom[i + 1, 1]
    }
    colour_index <- 1 + chrom[i, 2]
    colour_to_plot <- color_palette[colour_index]

    graphics::rect(xleft = xleft,
         xright = xrght,
         ybottom = 0,
         ytop = 1,
         col = colour_to_plot,
         border = NA)
  }
}
