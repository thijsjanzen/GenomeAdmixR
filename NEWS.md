# GenomeAdmixR 2.1.14

## Bug fixes

* `simulate_admixture()` ignored `ancestry_module(input_population = ...)` for
  single-population simulations and always started from a new random
  population.
* Selection in ancestry simulations ignored selected sites to the right of an
  individual's last junction, so a chromosome without crossovers never
  contributed to fitness. Several selected sites within one ancestry block
  were also undercounted.
* Selection site locations are now interpreted in Morgan, as documented;
  previously they were treated as relative positions when `morgan != 1`.
* `track_junctions = TRUE` returned values around 1.8e19 when a chromosome had
  no crossovers (unsigned underflow), and junction counts were off by one.
  `print.individual()` also reported one junction too few.
* Two-population ancestry simulations never filled in `junctions`, and did
  not scale tracked markers by `morgan`.
* Genotypes read from VCF files were assigned to the wrong rows, scrambling
  haplotypes and replacing about half of them with missing data.
* In sequence simulations without a `recombination_rate`, crossovers were
  placed proportional to marker position instead of marker spacing, so they
  were concentrated towards the end of the chromosome.
* Running sequence simulations with `num_threads > 1` could reuse the
  recombination map of an earlier simulation in the same R session, giving
  wrong results or crashing R.
* `calculate_allele_frequencies()` measured frequencies at rescaled positions
  when the last location was not 1.
* A single tracked marker was not scaled by `morgan`.
* When a population became fixed before `total_runtime`, the frequency table
  contained rows of zeros, the reported final generation was wrong, and the
  last simulated generation was discarded. Sequence simulations with mutation
  no longer stop early on fixation.
* Selecting more than one site in sequence data failed, and sites were
  checked against the tracked markers instead of the markers in the data.
* Multiplicative selection in sequence data gave fitness 0 when the selection
  matrix contained ignored rows.
* In two-population sequence simulations, the resampled starting population 1
  was not used, and the normalised substitution matrix was discarded.
* `calculate_ld()` failed when the population was smaller than
  `sampled_individuals`, and did not handle missing `markers`.
* Out-of-range allele values and missing markers in sequence data now give an
  informative error instead of reading or writing out of bounds.
* Populations with fewer than two individuals (or fewer than two individuals
  with positive fitness) now give an error instead of hanging.

## Changes

* Simulations are now reproducible with `set.seed()`, and give the same result
  for any value of `num_threads`.
* The package no longer changes the global RcppParallel thread settings;
  `num_threads < 1` uses all available cores.
* `calculate_ld()` samples individuals at random. Without `markers` it uses 100
  random markers, or up to 100 sites for sequence data.
* Because of the fixes above, simulation results differ from earlier versions.
