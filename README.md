# GenomeAdmixR  <img src="pics/GenomeAdmixR_2.png" align="right" width="180" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/GenomeAdmixR)](https://cran.r-project.org/package=GenomeAdmixR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GenomeAdmixR)](https://cran.r-project.org/package=GenomeAdmixR)
[![](http://cranlogs.r-pkg.org/badges/GenomeAdmixR)](https://cran.r-project.org/package=GenomeAdmixR)
[![R-CMD-check](https://github.com/thijsjanzen/GenomeAdmixR/workflows/R-CMD-check/badge.svg)](https://github.com/thijsjanzen/GenomeAdmixR/actions)

Branch|[![AppVeyor logo](pics/AppVeyor.png)]|[![Codecov logo](pics/Codecov.png)]
---|---|---
master|[![Build status](https://ci.appveyor.com/api/projects/status/vrfuo3dednjl52tr?svg=true)](https://ci.appveyor.com/project/thijsjanzen/genomeadmixr)|[![codecov.io](https://codecov.io/gh/thijsjanzen/GenomeAdmixR/branch/master/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/GenomeAdmixR)

# What is GenomeAdmixR?
A package to perfrom individual-based simulations of genome evolution with ancestry.

# Demonstration GenomeAdmixR
Thijs Janzen gave a presentation demonstrating GenomeAdmixR (then named isoSIM) at the R User Group at the University of Groningen, Groningen, The Netherlands. You can watch his presentation [here](https://p2gonew.serverspace.rug.nl/p2g/player/player.aspx?id=cxbKvM)  (audio starts after 1 min)

# More information
More information can be found in the accompanying paper, currently available as pre-print here : https://doi.org/10.1101/2020.10.19.343491

# Version history
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
