//
//  simulate_ancestry.cpp
//
//  Entry points for simulations that track local ancestry.
//  Created by Thijs Janzen on 28/02/2018.
//

#include <vector>

#include "Fish.h"
#include "random_functions.h"
#include "helper_functions.h"
#include "util.h"
#include "models.h"
#include "simulation.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

namespace {

// Written by Claude (Anthropic).
AncestryModel create_ancestry_model(const NumericMatrix& select,
                                    double morgan,
                                    bool multiplicative_selection) {
  AncestryModel model;
  model.morgan = morgan;
  // selection sites are given in Morgan, chromosomes are stored on [0, 1]
  model.selection_sites = convert_select_from_r(select, morgan);
  model.multiplicative_selection = multiplicative_selection;
  return model;
}

// R passes c(-1e6, -1e6) when there is no input population
// Written by Claude (Anthropic).
bool has_input_population(const NumericVector& input_population) {
  return input_population[0] > -1e4;
}

// Each individual is the offspring of two founders, which are drawn with
// probabilities 'founder_frequencies'.
// Written by Claude (Anthropic).
std::vector<Fish> create_founder_population(size_t population_size,
                                            const NumericVector& founder_frequencies,
                                            double morgan,
                                            rnd_t& random_generator) {
  std::vector<Fish> population;
  population.reserve(population_size);
  for (size_t i = 0; i < population_size; ++i) {
    const int first_founder = draw_random_founder(founder_frequencies,
                                                  random_generator);
    const int second_founder = draw_random_founder(founder_frequencies,
                                                   random_generator);
    const Fish first_parent(first_founder);
    const Fish second_parent(second_founder);
    population.push_back(mate(first_parent, second_parent, morgan,
                              random_generator));
  }
  return population;
}

// Written by Claude (Anthropic).
void add_founder_labels(const std::vector<Fish>& population,
                        std::vector<int>& founder_labels) {
  for (const auto& individual : population) {
    update_founder_labels(individual.chromosome1, founder_labels);
    update_founder_labels(individual.chromosome2, founder_labels);
  }
}

// Written by Claude (Anthropic).
void use_numbered_founder_labels(size_t number_of_founders,
                                 std::vector<int>& founder_labels) {
  for (size_t founder = 0; founder < number_of_founders; ++founder) {
    founder_labels.push_back(static_cast<int>(founder));
  }
}

// markers are given in Morgan, chromosomes are stored on [0, 1]
// Written by Claude (Anthropic).
std::vector<double> markers_on_relative_scale(const NumericVector& markers,
                                              double morgan) {
  const NumericVector scaled_markers = scale_markers(markers, morgan);
  return std::vector<double>(scaled_markers.begin(), scaled_markers.end());
}

// Written by Claude (Anthropic).
SimulationSettings create_settings(const std::vector<size_t>& population_sizes,
                                   double migration_rate,
                                   int total_runtime,
                                   bool track_frequency,
                                   const std::vector<double>& markers,
                                   bool track_junctions,
                                   bool verbose,
                                   int num_threads) {
  SimulationSettings settings;
  settings.population_sizes = population_sizes;
  settings.migration_rate = migration_rate;
  settings.total_runtime = total_runtime;
  settings.track_frequency = track_frequency;
  settings.tracked_markers = markers;
  settings.track_junctions = track_junctions;
  settings.summary_markers = markers;
  settings.verbose = verbose;
  settings.num_threads = num_threads;
  return settings;
}

// Written by Claude (Anthropic).
void print_message(bool verbose, const char* message) {
  if (verbose) {
    Rcout << message;
    force_output();
  }
}

}  // namespace

// Rewritten by Claude (Anthropic); the original was written by Thijs Janzen.
// [[Rcpp::export]]
List simulate_cpp(Rcpp::NumericVector input_population,
                  NumericMatrix select,
                  size_t pop_size,
                  size_t number_of_founders,
                  Rcpp::NumericVector starting_proportions,
                  size_t total_runtime,
                  double morgan,
                  bool verbose,
                  bool track_frequency,
                  NumericVector track_markers,
                  bool track_junctions,
                  bool multiplicative_selection,
                  int num_threads) {
  try {
    rnd_t random_generator(seed_from_r());
    AncestryModel model = create_ancestry_model(select, morgan,
                                                multiplicative_selection);

    std::vector<Fish> starting_population;
    if (has_input_population(input_population)) {
      print_message(verbose, "found input population, converting\n");
      const std::vector<Fish> input =
        convert_NumericVector_to_fishVector(input_population);
      add_founder_labels(input, model.founder_labels);
      starting_population = resize_population(input, pop_size,
                                              random_generator);
    } else {
      print_message(verbose, "starting generating random starting population\n");
      starting_population = create_founder_population(pop_size,
                                                      starting_proportions,
                                                      morgan,
                                                      random_generator);
      use_numbered_founder_labels(number_of_founders, model.founder_labels);
    }

    const SimulationSettings settings =
      create_settings({pop_size}, 0.0, total_runtime, track_frequency,
                      markers_on_relative_scale(track_markers, morgan),
                      track_junctions, verbose, num_threads);

    const auto result = run_simulation(model, settings,
                                       {starting_population},
                                       random_generator);
    print_message(verbose, "done simulating\n");

    return List::create(
      Named("population") = model.to_list(result.populations[0]),
      Named("frequencies") = result.frequencies,
      Named("initial_frequencies") = result.initial_frequencies,
      Named("final_frequencies") = result.final_frequencies,
      Named("junctions") = result.junctions);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// Rewritten by Claude (Anthropic); the original was written by Thijs Janzen.
// [[Rcpp::export]]
List simulate_migration_cpp(NumericVector input_population_1,
                            NumericVector input_population_2,
                            NumericMatrix select,
                            NumericVector pop_size,
                            NumericMatrix starting_frequencies,
                            int total_runtime,
                            double morgan,
                            bool verbose,
                            bool track_frequency,
                            NumericVector track_markers,
                            bool track_junctions,
                            bool multiplicative_selection,
                            double migration_rate,
                            int num_threads) {
  try {
    rnd_t random_generator(seed_from_r());
    AncestryModel model = create_ancestry_model(select, morgan,
                                                multiplicative_selection);

    const std::vector<size_t> population_sizes = {
      static_cast<size_t>(pop_size[0]),
      static_cast<size_t>(pop_size[1])
    };

    std::vector< std::vector<Fish> > starting_populations(2);
    if (has_input_population(input_population_1)) {
      const std::vector<Fish> input_1 =
        convert_NumericVector_to_fishVector(input_population_1);
      const std::vector<Fish> input_2 =
        convert_NumericVector_to_fishVector(input_population_2);
      add_founder_labels(input_1, model.founder_labels);
      add_founder_labels(input_2, model.founder_labels);
      starting_populations[0] = resize_population(input_1, population_sizes[0],
                                                  random_generator);
      starting_populations[1] = resize_population(input_2, population_sizes[1],
                                                  random_generator);
    } else {
      for (size_t index = 0; index < 2; ++index) {
        const NumericVector founder_frequencies = starting_frequencies(index, _);
        starting_populations[index] = create_founder_population(
                                          population_sizes[index],
                                          founder_frequencies,
                                          morgan,
                                          random_generator);
      }
      use_numbered_founder_labels(starting_frequencies.ncol(),
                                  model.founder_labels);
    }

    const SimulationSettings settings =
      create_settings(population_sizes, migration_rate, total_runtime,
                      track_frequency,
                      markers_on_relative_scale(track_markers, morgan),
                      track_junctions, verbose, num_threads);

    print_message(verbose, "starting simulation\n");
    const auto result = run_simulation(model, settings, starting_populations,
                                       random_generator);
    print_message(verbose, "done simulating\n");

    return List::create(
      Named("population_1") = model.to_list(result.populations[0]),
      Named("population_2") = model.to_list(result.populations[1]),
      Named("frequencies") = result.frequencies,
      Named("initial_frequencies") = result.initial_frequencies,
      Named("final_frequencies") = result.final_frequencies,
      Named("junctions") = result.junctions);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
