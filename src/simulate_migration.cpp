//
//  selection.cpp
//
//
//  Copyright Thijs Janzen 2018
//
//
#include <array>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <cmath>

#include <vector>
#include <algorithm>

#include <chrono>

#include "Fish.h"
#include "random_functions.h"
#include "helper_functions.h"
#include "util.h"

#include <RcppParallel.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

Fish draw_parent(const std::vector< Fish>& pop_1,
                 const std::vector< Fish>& pop_2,
                 double migration_rate,
                 bool use_selection,
                 const std::vector< double >& fitness_source,
                 const std::vector< double >& fitness_migr,
                 double max_fitness_source,
                 double max_fitness_migr,
                 int &index,
                 rnd_t& rndgen) {

  Fish parent;

  if (rndgen.uniform() < migration_rate) {
    // migration
    if(use_selection) {
      index = draw_prop_fitness(fitness_migr, max_fitness_migr, rndgen);
    } else {
      index = rndgen.random_number( static_cast<int>(pop_2.size()) );
    }
    parent = pop_2[index];
    // to ensure different indices for pop_1 and pop_2
    index += pop_1.size();
  } else {
    if(use_selection)  {
      index = draw_prop_fitness(fitness_source, max_fitness_source, rndgen);
    } else {
      index = rndgen.random_number( static_cast<int>(pop_1.size()) );
    }
    parent = pop_1[index];
  }
  return(parent);
}



std::vector< Fish > next_pop_migr_threaded(const std::vector< Fish>& pop_1,
                                           const std::vector< Fish>& pop_2,
                                           size_t pop_size,
                                           const std::vector< double >& fitness_source,
                                           const std::vector< double >& fitness_migr,
                                           const double& max_fitness_source,
                                           const double& max_fitness_migr,
                                           bool use_selection,
                                           double migration_rate,
                                           double size_in_morgan,
                                           int num_threads) {

  std::vector<Fish> new_generation(pop_size);

  rnd_t rndgen;
  int num_seeds = num_threads * 2; // tbb might re-start threads due to the load-balancer
  if (num_threads == -1) {
    num_seeds = 20;
  }
  std::vector< int > seed_values(num_seeds);

  for (int i = 0; i < num_seeds; ++i) {
    seed_values[i] = rndgen.random_number(INT_MAX); // large value
  }

  int seed_index = 0;
  std::mutex mutex;

  set_num_threads();

  tbb::parallel_for(
    tbb::blocked_range<unsigned>(0, pop_size),
    [&](const tbb::blocked_range<unsigned>& r) {

      thread_local rnd_t rndgen2(seed_values[seed_index]);
      {
        std::lock_guard<std::mutex> _(mutex);
        seed_index++;
        if (seed_index > num_seeds) { // just in case.
          for (int i = 0; i < num_seeds; ++i) {
            seed_values[i] = rndgen.random_number(INT_MAX);
          }
          seed_index = 0;
        }
      }

      for (unsigned i = r.begin(); i < r.end(); ++i) {
        int index1 = -1;
        int index2 = -1;
        Fish parent1 = draw_parent(pop_1, pop_2, migration_rate,
                                   use_selection,
                                   fitness_source, fitness_migr,
                                   max_fitness_source, max_fitness_migr,
                                   index1, rndgen2);
        Fish parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                   use_selection,
                                   fitness_source, fitness_migr,
                                   max_fitness_source, max_fitness_migr,
                                   index2, rndgen2);
        while (index1 == index2) {
          parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                use_selection,
                                fitness_source, fitness_migr,
                                max_fitness_source, max_fitness_migr,
                                index2, rndgen2);
        }

        new_generation[i] = mate(parent1, parent2, size_in_morgan, rndgen2);
      }
    });
  return new_generation;
}

std::vector< Fish > next_pop_migr(const std::vector< Fish>& pop_1,
                                           const std::vector< Fish>& pop_2,
                                           size_t pop_size,
                                           const std::vector< double >& fitness_source,
                                           const std::vector< double >& fitness_migr,
                                           const double& max_fitness_source,
                                           const double& max_fitness_migr,
                                           bool use_selection,
                                           double migration_rate,
                                           double size_in_morgan,
                                           int num_threads) {

  if (num_threads > 1) {
    return next_pop_migr_threaded(pop_1,
                                  pop_2,
                                  pop_size,
                                  fitness_source,
                                  fitness_migr,
                                  max_fitness_source,
                                  max_fitness_migr,
                                  use_selection,
                                  migration_rate,
                                  size_in_morgan,
                                  num_threads);
  }


  std::vector<Fish> new_generation(pop_size);

  rnd_t rndgen2; // calls get_seed

  for (unsigned i = 0; i < pop_size; ++i) {
    int index1 = -1;
    int index2 = -1;
    Fish parent1 = draw_parent(pop_1, pop_2, migration_rate,
                               use_selection,
                               fitness_source, fitness_migr,
                               max_fitness_source, max_fitness_migr,
                               index1, rndgen2);
    Fish parent2 = draw_parent(pop_1, pop_2, migration_rate,
                               use_selection,
                               fitness_source, fitness_migr,
                               max_fitness_source, max_fitness_migr,
                               index2, rndgen2);
    while (index1 == index2) {
      parent2 = draw_parent(pop_1, pop_2, migration_rate,
                            use_selection,
                            fitness_source, fitness_migr,
                            max_fitness_source, max_fitness_migr,
                            index2, rndgen2);
    }

    new_generation[i] = mate(parent1, parent2, size_in_morgan, rndgen2);
  }
  return new_generation;
}

std::vector< std::vector< Fish > > simulate_two_populations(
    const std::vector< Fish>& source_pop_1,
    const std::vector< Fish>& source_pop_2,
    const NumericMatrix& select,
    const NumericVector& pop_size,
    int total_runtime,
    double morgan,
    bool verbose,
    arma::mat& frequencies,
    bool track_frequency,
    const NumericVector& track_markers,
    bool track_junctions,
    std::vector<double>& junctions,
    bool multiplicative_selection,
    int num_alleles,
    const std::vector<int>& founder_labels,
    double migration_rate,
    int num_threads,
    rnd_t& rndgen) {

  std::vector< std::array<double, 5> > select_matrix;
  for (size_t i = 0; i < select.nrow(); ++i) {
    std::array<double, 5> row_entry;
    for (size_t j = 0; j < select.ncol(); ++j) {
       row_entry[j] = select(i, j);
    }
    select_matrix.push_back(row_entry);
  }

  bool use_selection = false;
  if (select_matrix[0][0] >= 0) use_selection = true;

  std::vector<Fish> pop_1 = source_pop_1;
  std::vector<Fish> pop_2 = source_pop_2;

  std::vector<double> fitness_pop_1(pop_1.size(), 0.0);
  std::vector<double> fitness_pop_2(pop_2.size(), 0.0);

  if (use_selection) {

    for (auto j : select_matrix) {
      if (j[4] < 0) break; // these entries are only for tracking, not for selection calculations

      double local_max_fitness = 0.0;
      for (int i = 1; i < 4; ++i) {
        if (j[i] > local_max_fitness) {
          local_max_fitness = j[i];
        }
      }
    }

    for (size_t i = 0; i < pop_1.size(); ++i) {
      fitness_pop_1[i] = calculate_fitness(pop_1[i],
                                           select_matrix,
                                           multiplicative_selection);
    }

    for (size_t i = 0; i < pop_2.size(); ++i) {
      fitness_pop_2[i] = calculate_fitness(pop_2[i],
                                           select_matrix,
                                           multiplicative_selection);
    }
  }

  double max_fitness_pop_1 = *std::max_element(fitness_pop_1.begin(), fitness_pop_1.end());
  double max_fitness_pop_2 = *std::max_element(fitness_pop_2.begin(), fitness_pop_2.end());


  int updateFreq = total_runtime / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(verbose) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }
  R_FlushConsole();

  for (int t = 0; t < total_runtime; ++t) {
    // Rcout << t << "\n"; force_output();
    if(track_frequency) {
      arma::mat local_mat = update_all_frequencies_tibble_dual_pop (pop_1,
                                                                    pop_2,
                                                                    track_markers,
                                                                    founder_labels,
                                                                    t,
                                                                    morgan);
      int num_founder_labels = founder_labels.size();
      int num_markers = track_markers.size();
      int local_mat_size = num_founder_labels * num_markers * 2;

      int start_add_marker = local_mat_size * t;

      for (int j = 0; j < local_mat_size; ++j) {
        for (int k = 0; k < 5; ++k) {
          frequencies(start_add_marker + j, k)  = local_mat(j, k);
        }
      }
    }

    std::vector<double> new_fitness_pop_1(pop_1.size(), 0.0);
    std::vector<double> new_fitness_pop_2(pop_2.size(), 0.0);

    std::vector<Fish> new_generation_pop_1 = next_pop_migr(pop_1, // resident
                                                           pop_2, // migrants
                                                           pop_size[0],
                                                           fitness_pop_1,
                                                           fitness_pop_2,
                                                           max_fitness_pop_1,
                                                           max_fitness_pop_2,
                                                           use_selection,
                                                           migration_rate,
                                                           morgan,
                                                           num_threads);

    std::vector<Fish> new_generation_pop_2 = next_pop_migr(pop_2,  // resident
                                                           pop_1,  // migrants
                                                           pop_size[1],
                                                           fitness_pop_2,
                                                           fitness_pop_1,
                                                           max_fitness_pop_2,
                                                           max_fitness_pop_1,
                                                           use_selection,
                                                           migration_rate,
                                                           morgan,
                                                           num_threads);
    pop_1 = new_generation_pop_1;
    pop_2 = new_generation_pop_2;

    if (use_selection) {
      for (size_t i = 0; i < pop_1.size(); ++i) {
        fitness_pop_1[i] = calculate_fitness(pop_1[i], select_matrix, multiplicative_selection);
      }
      for (size_t i = 0; i < pop_2.size(); ++i) {
        fitness_pop_2[i] = calculate_fitness(pop_2[i], select_matrix, multiplicative_selection);
      }
    }

    max_fitness_pop_1 = *std::max_element(fitness_pop_1.begin(),
                                          fitness_pop_1.end());
    max_fitness_pop_2 = *std::max_element(fitness_pop_2.begin(),
                                          fitness_pop_2.end());

    if (t % updateFreq == 0 && verbose) {
      Rcout << "**";
    }

    if (t > 2 && is_fixed(pop_1) && is_fixed(pop_2)) {
      if (verbose) Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
      force_output();
      std::vector< std::vector < Fish > > output;
      output.push_back(pop_1);
      output.push_back(pop_2);
      return(output);
    }

    Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcout << "\n";
  std::vector< std::vector< Fish > > output;
  output.push_back(pop_1);
  output.push_back(pop_2);
  return(output);
}

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
  rnd_t rndgen;

  std::vector< Fish > Pop_1;
  std::vector< Fish > Pop_2;
  int number_of_alleles = -1;
  std::vector<int> founder_labels;

  if (input_population_1[0] > -1e4) {
    Pop_1 = convert_NumericVector_to_fishVector(input_population_1);
    Pop_2 = convert_NumericVector_to_fishVector(input_population_2);

    for (auto it = Pop_1.begin(); it != Pop_1.end(); ++it) {
      update_founder_labels((*it).chromosome1, founder_labels);
      update_founder_labels((*it).chromosome2, founder_labels);
    }

    for (auto it = Pop_2.begin(); it != Pop_2.end(); ++it) {
      update_founder_labels((*it).chromosome1, founder_labels);
      update_founder_labels((*it).chromosome2, founder_labels);
    }

    number_of_alleles = founder_labels.size();

    if (Pop_1.size() != pop_size[0]) {
      // the populations have to be populated from the parents!
      std::vector< Fish > Pop_1_new;
      for(int j = 0; j < pop_size[0]; ++j) {
        int index = rndgen.random_number(Pop_1.size());
        Pop_1_new.push_back(Pop_1[index]);
      }
      Pop_1 = Pop_1_new;
      Pop_1_new.clear(); // free up memory
    }

    if (Pop_2.size() != pop_size[1]) {
      std::vector< Fish > Pop_2_new;
      for (int j = 0; j < pop_size[1]; ++j) {
        int index = rndgen.random_number(Pop_2.size());
        Pop_2_new.push_back(Pop_2[index]);
      }
      Pop_2 = Pop_2_new;
      Pop_2_new.clear(); // free up memory
    }
  } else {

    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < pop_size[j]; ++i) {
        NumericVector focal_freqs = starting_frequencies(j, _);

        int founder_1 = draw_random_founder(focal_freqs, rndgen);
        int founder_2 = draw_random_founder(focal_freqs, rndgen);

        Fish p1 = Fish( founder_1 );
        Fish p2 = Fish( founder_2 );

        if(j == 0) Pop_1.push_back(mate(p1, p2, morgan, rndgen));
        if(j == 1) Pop_2.push_back(mate(p1, p2, morgan, rndgen));
      }
    }
    for (int i = 0; i < starting_frequencies.ncol(); ++i) {
      founder_labels.push_back(i);
    }
    number_of_alleles = founder_labels.size();
  }

  int number_of_markers = track_markers.size();
  // 5 columns: time, loc, anc, type, population
  arma::mat frequencies_table(number_of_markers * number_of_alleles * total_runtime * 2, 5);
  arma::mat initial_frequencies = update_all_frequencies_tibble_dual_pop(Pop_1,
                                                                         Pop_2,
                                                                         track_markers,
                                                                         founder_labels,
                                                                         0,
                                                                         morgan);

  std::vector<double> junctions;
  std::vector< std::vector< Fish> > output_populations;
  if (verbose) {Rcout << "starting simulation\n"; force_output();}
  output_populations = simulate_two_populations(Pop_1,
                                                Pop_2,
                                                select,
                                                pop_size,
                                                total_runtime,
                                                morgan,
                                                verbose,
                                                frequencies_table,
                                                track_frequency,
                                                track_markers,
                                                track_junctions,
                                                junctions,
                                                multiplicative_selection,
                                                number_of_alleles,
                                                founder_labels,
                                                migration_rate,
                                                num_threads,
                                                rndgen);

  if (verbose) {Rcout << "done simulating\n"; force_output();}

  arma::mat final_frequencies = update_all_frequencies_tibble_dual_pop(output_populations[0],
                                                                       output_populations[1],
                                                                                         track_markers,
                                                                                         founder_labels,
                                                                                         total_runtime,
                                                                                         morgan);
  return List::create( Named("population_1") = convert_to_list(output_populations[0]),
                       Named("population_2") = convert_to_list(output_populations[1]),
                       Named("frequencies") = frequencies_table,
                       Named("initial_frequencies") = initial_frequencies,
                       Named("final_frequencies") = final_frequencies,
                       Named("junctions") = junctions);
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}