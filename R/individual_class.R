print.individual <- function(x, ...) {
  print("Individual with two Chromosomes")
  v1 <- paste("Chromosome 1:",
              length(x$chromosome1[, 1]) - 2,
              "junctions")
  v2 <- paste("Chromosome 2:",
              length(x$chromosome2[, 1]) - 2,
              "junctions")
  print(v1)
  print(v2)
}

print.population <- function(x, ...) {
  v1 <- paste("Population with",
              length(x),
              "individuals")
  print(v1)
}

plot.individual <- function(x, ...) {
  alleles_chrom1 <- unique(x$chromosome1[, 2])
  alleles_chrom2 <- unique(x$chromosome2[, 2])
  num_colors <- 1 + max(alleles_chrom1, alleles_chrom2)
  if(num_colors > 20) num_colors <- 20
  color_palette <- grDevices::rainbow(num_colors)

  par(mfrow = c(2, 1))
  par(mar = c(2, 2, 2, 2))
  plot(NA,
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")

  for (i in seq_along(x$chromosome1[, 1])) {
    xleft <- x$chromosome1[i, 1]
    xrght <- 1
    if (i < length(x$chromosome1[, 1])) {
      xrght <- x$chromosome1[i + 1, 1]
    }
    colour_index <- 1 + x$chromosome1[i, 2]
    colour_to_plot <- color_palette[colour_index]

    rect(xleft = xleft,
         xright = xrght,
         ybottom = 0,
         ytop = 1,
         col = colour_to_plot,
         border = NA)
  }

  plot(NA,
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")

  for (i in seq_along(x$chromosome2[, 1])) {
    xleft <- x$chromosome2[i, 1]
    xrght <- 1
    if (i < length(x$chromosome2[, 1])) {
      xrght <- x$chromosome2[i + 1, 1]
    }
    colour_index <- 1 + x$chromosome2[i, 2]
    colour_to_plot <- color_palette[colour_index]

    rect(xleft = xleft,
         xright = xrght,
         ybottom = 0,
         ytop = 1,
         col = colour_to_plot,
         border = NA)
  }
}

plot_chromosome <- function(chrom, xmin = 0, xmax = 1) {
  alleles <- unique(chrom[, 2])
  num_colors <- 1 + max(alleles)
  if(num_colors > 20) num_colors <- 20
  color_palette <- grDevices::rainbow(num_colors)

  plot(NA,
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

    rect(xleft = xleft,
         xright = xrght,
         ybottom = 0,
         ytop = 1,
         col = colour_to_plot,
         border = NA)
  }
}

create_pop_class <- function(pop) {

  set_indiv_class <- function(indiv) {
    class(indiv) <- "individual"
    indiv
  }
  pop <- lapply(pop, set_indiv_class)
  class(pop) <- "population"
  return(pop)
}

verify_individual <- function(indiv) {

  if(!is(indiv, "individual")) return(FALSE)

  if(indiv$chromosome1[1,1] != 0) {
     cat("Chromosome doesn't start at 0\n")
     return(FALSE)
  }
  if(tail(indiv$chromosome1,1)[2] != -1) {
    cat("Chromosome doesn't end with -1\n")
    return(FALSE)
  }

  if(max(abs(indiv$chromosome1[,2])) > 10000) {
    cat("Memory error recorded in chromosome\n")
    return(FALSE)
  }

  if(indiv$chromosome2[1,1] != 0) {
    cat("Chromosome doesn't start at 0\n")
    return(FALSE)
  }

  if(tail(indiv$chromosome2,1)[2] != -1) {
    cat("Chromosome doesn't end with -1\n")
    return(FALSE)
  }

  if(max(abs(indiv$chromosome2[,2])) > 10000) {
    cat("Memory error recorded in chromosome\n")
    return(FALSE)
  }

  return(TRUE)
}

verify_population <- function(pop) {

  if(!is(pop, "population"))  {
    if(!is(pop$population, "population")) {
      return(FALSE)
    } else {
      return(verify_population(pop$population))
    }
  }
  v <- unlist(lapply(pop, verify_individual))
  if(sum(v) != length(v)) return(FALSE)

  return(TRUE)
}


findtype <- function(chrom, pos) {

  b <- which(chrom[, 1] > pos)
  chromtype <- chrom[b[1] - 1, 2]

  return(chromtype[[1]])
}