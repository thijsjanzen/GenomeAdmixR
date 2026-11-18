//
//  simulate_sequence.cpp
//
//  Entry points for simulations using sequence data.
//  Created by Thijs Janzen on 28/02/2018.
//

#include <numeric>
#include <vector>

#include "Fish_emp.h"
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
bool has_user_recombination_map(const NumericVector& recombination_map,
                                size_t number_of_locations) {
  return recombination_map.size() > 0 &&
         static_cast<size_t>(recombination_map.size()) == number_of_locations;
}

// A user supplied map also sets the length of the chromosome (its sum, in
// Morgan). Otherwise crossovers are proportional to the distance between
// adjacent locations.
// Written by Claude (Anthropic).
void set_recombination_map(SequenceModel& model,
                           const NumericVector& recombination_map) {
  if (has_user_recombination_map(recombination_map, model.locations.size())) {
    const std::vector<double> map_values(recombination_map.begin(),
                                         recombination_map.end());
    model.recombination_map = emp_genome(map_values);
    model.morgan = std::accumulate(map_values.begin(), map_values.end(), 0.0);
  } else {
    model.recombination_map = emp_genome::from_positions(model.locations);
  }
}

// Written by Claude (Anthropic).
std::vector< std::vector<double> > to_substitution_matrix(
    const NumericMatrix& substitution_matrix) {
  std::vector< std::vector<double> > output(4, std::vector<double>(4, 0.0));
  for (int from_base = 0; from_base < 4; ++from_base) {
    for (int to_base = 0; to_base < 4; ++to_base) {
      output[from_base][to_base] = substitution_matrix(from_base, to_base);
    }
  }
  return output;
}

// Written by Claude (Anthropic).
void set_mutation(SequenceModel& model,
                  double mutation_rate,
                  const NumericMatrix& substitution_matrix) {
  model.mutation_rate = mutation_rate;
  if (mutation_rate > 0) {
    model.substitution_matrix = to_substitution_matrix(substitution_matrix);
  } else {
    model.substitution_matrix.assign(4, std::vector<double>(4, 0.0));
  }
}

// Written by Claude (Anthropic).
SequenceModel create_sequence_model(const NumericVector& locations,
                                    const NumericMatrix& select,
                                    double morgan,
                                    bool multiplicative_selection,
                                    double mutation_rate,
                                    const NumericMatrix& substitution_matrix,
                                    const NumericVector& recombination_map) {
  SequenceModel model;
  model.locations.assign(locations.begin(), locations.end());
  model.selection_sites = convert_select_from_r(select);
  model.multiplicative_selection = multiplicative_selection;
  model.morgan = morgan;
  set_recombination_map(model, recombination_map);
  set_mutation(model, mutation_rate, substitution_matrix);
  return model;
}

// R passes a matrix filled with -1e6 when there is no input population
// Written by Claude (Anthropic).
bool has_input_population(const NumericMatrix& input_population) {
  return input_population(0, 0) > -1e4;
}

// Written by Claude (Anthropic).
std::vector<Fish_emp> create_starting_population(const NumericMatrix& input_population,
                                                 size_t population_size,
                                                 rnd_t& random_generator) {
  const std::vector<Fish_emp> input =
    convert_numeric_matrix_to_fish_vector(input_population);
  return resize_population(input, population_size, random_generator);
}

// Written by Claude (Anthropic).
SimulationSettings create_settings(const std::vector<size_t>& population_sizes,
                                   double migration_rate,
                                   int total_runtime,
                                   bool track_frequency,
                                   const std::vector<double>& tracked_markers,
                                   const std::vector<double>& summary_markers,
                                   bool verbose,
                                   int num_threads) {
  SimulationSettings settings;
  settings.population_sizes = population_sizes;
  settings.migration_rate = migration_rate;
  settings.total_runtime = total_runtime;
  settings.track_frequency = track_frequency;
  settings.tracked_markers = tracked_markers;
  settings.track_junctions = false;  // sequence data has no junctions
  settings.summary_markers = summary_markers;
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
List simulate_emp_cpp(const Rcpp::NumericMatrix& input_population,
                      const Rcpp::NumericVector& marker_positions_R,
                      const Rcpp::NumericMatrix& select,
                      size_t pop_size,
                      size_t total_runtime,
                      double morgan,
                      bool verbose,
                      bool track_frequency,
                      const Rcpp::NumericVector& track_markers_R,
                      bool multiplicative_selection,
                      double mutation_rate,
                      const Rcpp::NumericMatrix& substitution_matrix_R,
                      int num_threads,
                      const Rcpp::NumericVector& recombination_map) {
  try {
    rnd_t random_generator(seed_from_r());
    const SequenceModel model = create_sequence_model(marker_positions_R,
                                                      select, morgan,
                                                      multiplicative_selection,
                                                      mutation_rate,
                                                      substitution_matrix_R,
                                                      recombination_map);

    if (!has_input_population(input_population)) {
      Rcpp::stop("can not run without input data");
    }
    const std::vector<Fish_emp> starting_population =
      create_starting_population(input_population, pop_size, random_generator);

    // initial and final frequencies are reported for all locations
    const std::vector<double> tracked_markers(track_markers_R.begin(),
                                              track_markers_R.end());
    const SimulationSettings settings =
      create_settings({pop_size}, 0.0, total_runtime, track_frequency,
                      tracked_markers, model.locations, verbose, num_threads);

    const auto result = run_simulation(model, settings,
                                       {starting_population},
                                       random_generator);

    return List::create(
      Named("population") = model.to_list(result.populations[0]),
      Named("frequencies") = result.frequencies,
      Named("initial_frequencies") = result.initial_frequencies,
      Named("final_frequencies") = result.final_frequencies);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// Rewritten by Claude (Anthropic); the original was written by Thijs Janzen.
// [[Rcpp::export]]
List simulate_migration_emp_cpp(const NumericMatrix& input_population_1,
                                const NumericMatrix& input_population_2,
                                const NumericVector& marker_positions_R,
                                NumericMatrix select,
                                const NumericVector& pop_sizes,
                                int total_runtime,
                                double morgan,
                                bool verbose,
                                bool track_frequency,
                                const NumericVector& track_markers_R,
                                bool multiplicative_selection,
                                double migration_rate,
                                double mutation_rate,
                                const NumericMatrix& substitution_matrix_R,
                                int num_threads,
                                const NumericVector& recombination_map) {
  try {
    rnd_t random_generator(seed_from_r());
    const SequenceModel model = create_sequence_model(marker_positions_R,
                                                      select, morgan,
                                                      multiplicative_selection,
                                                      mutation_rate,
                                                      substitution_matrix_R,
                                                      recombination_map);

    if (!has_input_population(input_population_1)) {
      Rcpp::stop("empty populations, sequence simulations need input data");
    }
    const std::vector<size_t> population_sizes = {
      static_cast<size_t>(pop_sizes[0]),
      static_cast<size_t>(pop_sizes[1])
    };
    std::vector< std::vector<Fish_emp> > starting_populations(2);
    starting_populations[0] = create_starting_population(input_population_1,
                                                         population_sizes[0],
                                                         random_generator);
    starting_populations[1] = create_starting_population(input_population_2,
                                                         population_sizes[1],
                                                         random_generator);

    const std::vector<double> tracked_markers(track_markers_R.begin(),
                                              track_markers_R.end());
    const SimulationSettings settings =
      create_settings(population_sizes, migration_rate, total_runtime,
                      track_frequency, tracked_markers, tracked_markers,
                      verbose, num_threads);

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
      Named("junctions") = std::vector<double>());
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
