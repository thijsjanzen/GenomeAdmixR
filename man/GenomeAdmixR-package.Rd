\name{GenomeAdmixR-package}
\alias{GenomeAdmixR-package}
\alias{GenomeAdmixR}
\docType{package}
\title{
 Simulate Admixture of Genomes
}
\description{
   Individual-based simulations forward in time, simulating how
patterns in ancestry along the genome change after admixture. The simulation
assumes Wright-Fisher dynamics, e.g. random mating and non-overlapping generations. In the simulation, instead of specific alleles, local ancestry is tracked, thus assuming that local molecular data can always be uniquely traced back to one of the founding individuals (populations). The package provides functionality to perform such simulations, but also to perform post-hoc statistical analyses and to visualize the obtained results.

   Version 2.0.1 - Moved migration outside the modules \cr
   Version 2.0  - Added ancestry_module and sequence_module to distinguish between
                  implementations of the model \cr
   Version 1.2  - Added example sequencing data \cr
   Version 1.2  - Added the option to load sequence data for admixing \cr
   Version 1.1  - Fixed a minor bug with plot_joyplot_frequencies \cr
   Version 1.1  - Improved tests \cr
   Version 1.1  - Improved recombination code (again) \cr
   Version 1.0  - Release associated with bioRxiv submission, to be found here: https://doi.org/10.1101/2020.10.19.343491 \cr
   Version 0.66 - Improved recombination code, about twice as fast \cr
   Version 0.65 - Added testing and added logo \cr
   Version 0.64 - Reduced cyclomatic complexity \cr
   Version 0.63 - Updated random number generation \cr
   Version 0.62 - Updated to Roxygen \cr
   Version 0.61 - Added plot_over_time \cr
   Version 0.60 - Added admixture with migration \cr
   Version 0.59 - Updated frequency code under the hood \cr
   Version 0.58 - Renamed to GenomeAdmixR \cr
   Version 0.58 - Collapsed and improved many functions \cr
   Version 0.57 - Added function to generate admixed individuals \cr
   Version 0.56 - Added starting frequencies to 'simulate_admixture' \cr
   Version 0.55 - Extended 'calculate_marker_frequency' to handle a vector of locations \cr
   Version 0.55 - Increased accuracy of choosing a random position for recombination, this should prevent the rare bug fixed in version 0.54 \cr
   Version 0.54 - Fixed a MAJOR bug regarding recombination: in rare cases, a crossover position could be picked on an existing junction, due to the limited number of digits in uniform() \cr
   Version 0.54 - Improved plot_difference_frequencies to handle modified input \cr
   Version 0.53 - Added multiplicative_selection \cr
   Version 0.52 - Added plot_difference_frequencies \cr
   Version 0.51 - Added tajima's d calculation \cr
   Version 0.50 - Added simulated_admixture until \cr
   Version 0.49 - Added 'simulate' to cpp \cr
   Version 0.48 - Added a general 'simulate' function \cr
   Version 0.47 - Changed the effect of migration \cr
   Version 0.46 - Added joyplot & increase_ancestor \cr
   Version 0.45 - Removed create_two_populations \cr
   Version 0.44 - Added tracking regions \cr
   Version 0.43 - Fixed bugs in select_population \cr
   Version 0.42 - Added initial and final frequency tables \cr
   Version 0.41 - Added multiple marker support \cr
   Version 0.40 - Collapsed selection functions \cr
   Version 0.39 - Added support for non-additive selection \cr
   Version 0.38 - Added track frequencies \cr
   Version 0.37 - Removed selection on regions \cr
   Version 0.36 - Added progress_bar option \cr
   Version 0.35 - Added calculate_marker_frequency \cr
   Version 0.34 - Added selection_markers \cr
   Version 0.33 - Fixed bugs in selection \cr
   Version 0.32 - Moved Fish.h code to Fish.cpp \cr
   Version 0.31 - Changed random number generator to R based \cr
   Version 0.30 - Added Recombination = 1 code \cr
   Version 0.29 - Changed internal junction representation: removed .left \cr
   Version 0.28 - Reverted to Agner Fog Random number generation \cr
   Version 0.27 - Speed up return types \cr
   Version 0.26 - Added class verification code \cr
   Version 0.25 - Squashed plotting bug \cr
   Version 0.24 - Removed Output.cpp \cr
   Version 0.23 - Removed number_of_founders from calc_allele_spectrum \cr
   Version 0.22 - Added save and load functions \cr
   Version 0.21 - Changed random-seed management \cr
   Version 0.20 - Removed superfluous code \cr
   Version 0.19 - Removed number_of_founders from Fst and LD code \cr
   Version 0.18 - Start of tracking changes \cr
}
\author{
Thijs Janzen
Maintainer: (thijsjanzen@gmail.com)
}
\references{
Janzen, T.,  Diaz, F. (2020) Individual-based simulations of genome evolution with ancestry: the GenomeAdmixR R package. bioRxiv 2020.10.19.343491; doi: https://doi.org/10.1101/2020.10.19.343491
}
