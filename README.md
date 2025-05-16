# GenomeAdmixR  <img src="pics/GenomeAdmixR_2.png" align="right" width="180" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/GenomeAdmixR)](https://cran.r-project.org/package=GenomeAdmixR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GenomeAdmixR)](https://cran.r-project.org/package=GenomeAdmixR)
[![](http://cranlogs.r-pkg.org/badges/GenomeAdmixR)](https://cran.r-project.org/package=GenomeAdmixR)
[![R-CMD-check](https://github.com/thijsjanzen/GenomeAdmixR/workflows/R-CMD-check/badge.svg)](https://github.com/thijsjanzen/GenomeAdmixR/actions)

Branch|[![AppVeyor logo](pics/AppVeyor.png)]|[![Codecov logo](pics/Codecov.png)]
---|---|---
master|[![Build status](https://ci.appveyor.com/api/projects/status/vrfuo3dednjl52tr?svg=true)](https://ci.appveyor.com/project/thijsjanzen/genomeadmixr)|[![codecov.io](https://codecov.io/gh/thijsjanzen/GenomeAdmixR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/thijsjanzen/GenomeAdmixR/)

# What is GenomeAdmixR?
A package to perfrom individual-based simulations of genome evolution with ancestry.

# More information
More information can be found in the accompanying paper :  https://doi.org/10.1111/2041-210X.13612

# Version history
Version 2.1.10 - Fixed memory bug, improved documentation <br />
Version 2.1.9  - updated tbb::task_scheduler_init to tbb::global_control <br />
Version 2.1.7  - Improve documentation <br />
Version 2.1.6  - check classes with inherits <br />
Version 2.1.5  - Removed debugging output <br />
Version 2.1.4  - Only output when verbose = TRUE <br />
Version 2.1.3  - Changed DOI link in description <br />
Version 2.1.2  - Improved testing <br />
Version 2.1.1  - Removed GNU make dependency <br />
Version 2.1    - Removed error in calculate_allele_frequency <br />
Version 2.0.1 - Moved migration outside the modules <br />
Version 2.0  - Added ancestry_module and sequence_module to distinguish between implementations of the model <br />
Version 1.2  - Added example sequencing data <br />
Version 1.2  - Added the option to load sequence data for admixing <br />
Version 1.1  - Fixed a minor bug with plot_joyplot_frequencies <br />
Version 1.1  - Improved tests <br />
Version 1.1  - Improved recombination code (again) <br />
Version 1.00 - Release associated with bioRxiv submission, to be found here: https://doi.org/10.1101/2020.10.19.343491 <br />
Version 0.66 - Improved recombination code, about twice as fast <br />
Version 0.65 - Added testing and added logo <br />
Version 0.64 - Reduced cyclomatic complexity <br />
Version 0.63 - Updated random number generation <br />
Version 0.62 - Updated to Roxygen <br />
Version 0.61 - Added plot_over_time <br />
Version 0.60 - Added admixture with migration <br />
Version 0.59 - Updated underlying code tracking frequencies <br />
Version 0.58 - Removed many old functions, and improved usability for many existing functions <br />
Version 0.58 - Renamed to GenomeAdmixR <br />
Version 0.57 - Added function to generate admixed individuals <br />
Version 0.56 - Added starting frequencies to 'simulate_admixture' <br />
Version 0.55 - extended 'calculate_marker_frequency' to handle a vector of locations <br />
Version 0.55 - increased accuracy of choosing a random position for recombination, this should prevent the rare bug fixed in version 0.54 <br />
Version 0.54 - Fixed a MAJOR bug regarding recombination: in rare cases, a crossover position could be picked on an existing junction, due to the limited number of digits in uniform() <br />
Version 0.54 - Improved plot_difference_frequencies to handle modified input <br />
Version 0.53 - Added multiplicative_selection <br />
Version 0.52 - Added plot_difference_frequencies <br />
Version 0.51 - Added tajima's d calculation <br />
Version 0.50 - added simulated_admixture until <br />
Version 0.49 - Added 'simulate' to cpp <br />
Version 0.48 - Added a general 'simulate' function <br />
Version 0.47 - Changed the effect of migration <br />
Version 0.46 - Added joyplot & increase_ancestor <br />
Version 0.45 - Removed create_two_populations <br />
Version 0.44 - Added tracking regions <br />
Version 0.43 - Fixed bugs in select_population <br />
Version 0.42 - Added initial and final frequency tables <br />
Version 0.41 - Added multiple marker support <br />
Version 0.40 - Collapsed selection functions <br />
Version 0.39 - Added support for non-additive selection <br />
Version 0.38 - Added track frequencies <br />
Version 0.37 - Removed selection on regions <br />
Version 0.36 - Added progress_bar option <br />
Version 0.35 - Added calculate_marker_frequency <br />
Version 0.34 - Added selection_markers <br />
Version 0.33 - Fixed bugs in selection <br />
Version 0.32 - Moved Fish.h code to Fish.cpp <br />
Version 0.31 - Changed random number generator to R based <br />
Version 0.30 - Added Recombination = 1 code <br />
Version 0.29 - Changed internal junction representation: removed .left <br />
Version 0.28 - Reverted to Agner Fog Random number generation <br />
Version 0.27 - Speed up return types <br />
Version 0.26 - Added class verification code <br />
Version 0.25 - Squashed plotting bug <br />
Version 0.24 - Removed Output.cpp <br />
Version 0.23 - Removed number_of_founders from calc_allele_spectrum <br />
Version 0.22 - Added save and load functions <br />
Version 0.21 - Changed random-seed management <br />
Version 0.20 - Removed superfluous code <br />
Version 0.19 - Removed number_of_founders from Fst and LD code <br />
Version 0.18 - Start of tracking changes <br />
