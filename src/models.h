//
//  models.h
//
//  The two ways GenomeAdmixR stores genomes, written as Models for the generic
//  simulation in simulation.h.
//

#pragma once

#include <array>
#include <utility>
#include <vector>

#include "Fish.h"
#include "Fish_emp.h"
#include "random_functions.h"
#include "helper_functions.h"


// Local ancestry along a chromosome of 'morgan' Morgan. Chromosomes are stored
// as junctions on [0, 1], so marker and selection locations are relative
// positions in [0, 1].
struct AncestryModel {
  using individual_type = Fish;

  double morgan = 1.0;
  std::vector<int> founder_labels;
  std::vector< std::array<double, 5> > selection_sites;
  bool multiplicative_selection = true;

  individual_type reproduce(const individual_type& first_parent,
                            const individual_type& second_parent,
                            rnd_t& random_generator) const {
    return mate(first_parent, second_parent, morgan, random_generator);
  }

  bool use_selection() const {
    return !selection_sites.empty();
  }

  double fitness(const individual_type& individual) const {
    return calculate_fitness(individual, selection_sites,
                             multiplicative_selection);
  }

  int number_of_alleles() const {
    return static_cast<int>(founder_labels.size());
  }

  arma::mat allele_frequencies(const std::vector<individual_type>& population,
                               const std::vector<double>& markers,
                               int generation) const {
    return update_all_frequencies_tibble(population, markers, founder_labels,
                                         generation, morgan);
  }

  double count_junctions(const std::vector<individual_type>& population) const {
    return total_junctions(population);
  }

  bool may_stop_when_fixed() const {
    return true;
  }

  Rcpp::List to_list(const std::vector<individual_type>& population) const {
    return convert_to_list(population);
  }
};


// Sequence data: a base (0 = missing, 1 to 4 = a, c, t, g) at every location.
// Marker and selection locations are locations present in the data.
struct SequenceModel {
  using individual_type = Fish_emp;

  double morgan = 1.0;
  std::vector<double> locations;
  emp_genome recombination_map;
  std::vector< std::array<double, 5> > selection_sites;
  bool multiplicative_selection = true;
  double mutation_rate = 0.0;
  std::vector< std::vector<double> > substitution_matrix;

  individual_type reproduce(const individual_type& first_parent,
                            const individual_type& second_parent,
                            rnd_t& random_generator) const {
    // The gametes are drawn in two separate statements, so the order in which
    // random numbers are used does not depend on the compiler. The second
    // parent's gamete is drawn first, which matches earlier results on GCC.
    auto gamete_of_second_parent = second_parent.gamete(morgan,
                                                        random_generator,
                                                        recombination_map);
    auto gamete_of_first_parent = first_parent.gamete(morgan,
                                                      random_generator,
                                                      recombination_map);
    individual_type offspring(std::move(gamete_of_first_parent),
                              std::move(gamete_of_second_parent));
    if (mutation_rate > 0) {
      mutate(offspring, substitution_matrix, mutation_rate, random_generator);
    }
    return offspring;
  }

  bool use_selection() const {
    return !selection_sites.empty();
  }

  double fitness(const individual_type& individual) const {
    return calculate_fitness(individual, selection_sites, locations,
                             multiplicative_selection);
  }

  int number_of_alleles() const {
    return 5;  // 0 (missing), a, c, t, g
  }

  arma::mat allele_frequencies(const std::vector<individual_type>& population,
                               const std::vector<double>& markers,
                               int generation) const {
    return update_all_frequencies_tibble(population, markers, locations,
                                         generation, morgan);
  }

  double count_junctions(const std::vector<individual_type>&) const {
    return 0.0;  // sequence data has no junctions
  }

  // with mutation, a fixed population can still change
  bool may_stop_when_fixed() const {
    return mutation_rate <= 0;
  }

  Rcpp::List to_list(const std::vector<individual_type>& population) const {
    return convert_to_list(population, locations);
  }
};
