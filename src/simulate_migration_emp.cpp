//
//  selection.cpp
//
//
//  Copyright Thijs Janzen 2018
//
//
#include <vector>
#include <cstdlib>
#include <numeric>
#include <cmath>

#include <vector>
#include <algorithm>

#include "Fish_emp.h"
#include "random_functions.h"
#include "helper_functions.h"
#include "util.h"

#include <RcppParallel.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

Fish_emp draw_parent(const std::vector< Fish_emp >& pop_1,
                     const std::vector< Fish_emp >& pop_2,
                     double migration_rate,
                     bool use_selection,
                     const std::vector< double >& fitness_source,
                     const std::vector< double >& fitness_migr,
                     double max_fitness_source,
                     double max_fitness_migr,
                     int &index,
                     rnd_t& rndgen) {

  Fish_emp parent;

  if (rndgen.uniform() < migration_rate) {
    // migration
    if(use_selection) {
      index = draw_prop_fitness(fitness_migr, max_fitness_migr, rndgen);
    } else {
      index = rndgen.random_number( (int)pop_2.size() );
    }
    parent = pop_2[index];
    index = index + pop_1.size();
    // to ensure different indices for pop_1 and pop_2
  } else {
    if(use_selection)  {
      index = draw_prop_fitness(fitness_source, max_fitness_source, rndgen);
    } else {
      index = rndgen.random_number( (int)pop_1.size() );
    }
    parent = pop_1[index];
  }
  return(parent);
}



std::vector< Fish_emp > next_pop_migr(const std::vector< Fish_emp >& pop_1,
                                      const std::vector< Fish_emp >& pop_2,
                                      const std::vector< double >& marker_positions,
                                      size_t pop_size,
                                      std::vector< double > fitness_source,
                                      std::vector< double > fitness_migr,
                                      double max_fitness_source,
                                      double max_fitness_migr,
                                      const NumericMatrix& select,
                                      bool use_selection,
                                      bool multiplicative_selection,
                                      double migration_rate,
                                      double size_in_morgan,
                                      double mutation_rate,
                                      const std::vector<std::vector<double>>& substitution_matrix,
                                      const emp_genome& emp_gen,
                                      int num_threads) {

  std::vector<Fish_emp> new_generation(pop_size);

  if (num_threads == 1) {
    rnd_t rndgen2;

    for (size_t i = 0; i < pop_size; ++i) {

      int index1, index2;
      Fish_emp parent1 = draw_parent(pop_1, pop_2, migration_rate,
                                     use_selection,
                                     fitness_source, fitness_migr,
                                     max_fitness_source, max_fitness_migr,
                                     index1,
                                     rndgen2);
      Fish_emp parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                     use_selection,
                                     fitness_source, fitness_migr,
                                     max_fitness_source, max_fitness_migr,
                                     index2,
                                     rndgen2);
      while (index1 == index2) {
        parent2 = draw_parent(pop_1, pop_2, migration_rate,
                              use_selection,
                              fitness_source, fitness_migr,
                              max_fitness_source, max_fitness_migr,
                              index2,
                              rndgen2);
      }

      new_generation[i] = Fish_emp(parent1.gamete(size_in_morgan, rndgen2, emp_gen),
                                   parent2.gamete(size_in_morgan, rndgen2, emp_gen));

      if (mutation_rate > 0)
        mutate(new_generation[i], substitution_matrix, mutation_rate, rndgen2);
    }
  } else {

    int seed_index = 0;
    std::mutex mutex;
    int num_seeds = num_threads * 10; // tbb might re-start threads due to the load-balancer
    if (num_threads == -1) {
      num_seeds = 200;
    }
    std::vector< int > seed_values(num_seeds);

    {
      rnd_t rndgen;
      for (int i = 0; i < num_seeds; ++i) {
        seed_values[i] = rndgen.random_number(INT_MAX); // large value
      }
    }

    set_num_threads();

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, pop_size),
      [&](const tbb::blocked_range<unsigned>& r) {

        thread_local rnd_t rndgen2(seed_values[seed_index]);
        {
          std::lock_guard<std::mutex> _(mutex);
          seed_index++;
          if (seed_index >= num_seeds) { // just in case.
            for (int i = 0; i < num_seeds; ++i) {
              seed_values[i] = rndgen2.random_number(INT_MAX);
            }
            seed_index = 0;
          }
        }

        for (unsigned i = r.begin(); i < r.end(); ++i) {

          int index1, index2;
          Fish_emp parent1 = draw_parent(pop_1, pop_2, migration_rate,
                                         use_selection,
                                         fitness_source, fitness_migr,
                                         max_fitness_source, max_fitness_migr,
                                         index1,
                                         rndgen2);
          Fish_emp parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                         use_selection,
                                         fitness_source, fitness_migr,
                                         max_fitness_source, max_fitness_migr,
                                         index2,
                                         rndgen2);
          while (index1 == index2) {
            parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                  use_selection,
                                  fitness_source, fitness_migr,
                                  max_fitness_source, max_fitness_migr,
                                  index2,
                                  rndgen2);
          }

          new_generation[i] = Fish_emp(parent1.gamete(size_in_morgan, rndgen2, emp_gen),
                                       parent2.gamete(size_in_morgan, rndgen2, emp_gen));

          if (mutation_rate > 0)
            mutate(new_generation[i], substitution_matrix, mutation_rate, rndgen2);
        }
      });
  }
  return new_generation;
}

std::vector< std::vector< Fish_emp > > simulate_two_populations(
    const std::vector< Fish_emp>& source_pop_1,
    const std::vector< Fish_emp>& source_pop_2,
    const std::vector<double>& marker_positions,
    const NumericMatrix& select_r,
    const std::vector<size_t>& pop_size,
    int total_runtime,
    double morgan,
    bool verbose,
    arma::mat& frequencies,
    bool track_frequency,
    const std::vector<double>& track_markers,
    bool multiplicative_selection,
    int num_alleles,
    const std::vector<int>& founder_labels,
    double migration_rate,
    double mutation_rate,
    const std::vector<std::vector<double>>& substitution_matrix,
    rnd_t& rndgen,
    const emp_genome& emp_gen,
    int num_threads) {

  bool use_selection = false;
  if (select_r(0, 0) >= 0) use_selection = true;

  std::vector<Fish_emp> pop_1 = source_pop_1;
  std::vector<Fish_emp> pop_2 = source_pop_2;

  if (pop_1.empty() || pop_2.empty()) {
    throw "empty populations";
  }

  std::vector<double> fitness_pop_1(pop_1.size(), 0.0);
  std::vector<double> fitness_pop_2(pop_2.size(), 0.0);

  std::vector< std::array< double, 5 >> select;

  if (use_selection) {
    select = convert_select_from_r(select_r);
    for (size_t i = 0; i < pop_1.size(); ++i) {
      fitness_pop_1[i] = calculate_fitness(pop_1[i], select,
                                           marker_positions,
                                           multiplicative_selection);
    }

    for (size_t i = 0; i < pop_2.size(); ++i) {
      fitness_pop_2[i] = calculate_fitness(pop_2[i],
                                           select,
                                           marker_positions,
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
                                                                    marker_positions,
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

    std::vector<Fish_emp> new_generation_pop_1 = next_pop_migr(pop_1, // resident
                                                               pop_2, // migrants
                                                               marker_positions,
                                                               pop_size[0],
                                                               fitness_pop_1,
                                                               fitness_pop_2,
                                                               max_fitness_pop_1,
                                                               max_fitness_pop_2,
                                                               select_r,
                                                               use_selection,
                                                               multiplicative_selection,
                                                               migration_rate,
                                                               morgan,
                                                               mutation_rate,
                                                               substitution_matrix,
                                                               emp_gen,
                                                               num_threads);

    std::vector<Fish_emp> new_generation_pop_2 = next_pop_migr(pop_2,  // resident
                                                               pop_1,  // migrants
                                                               marker_positions,
                                                               pop_size[1],
                                                               fitness_pop_2,
                                                               fitness_pop_1,
                                                               max_fitness_pop_2,
                                                               max_fitness_pop_1,
                                                               select_r,
                                                               use_selection,
                                                               multiplicative_selection,
                                                               migration_rate,
                                                               morgan,
                                                               mutation_rate,
                                                               substitution_matrix,
                                                               emp_gen,
                                                               num_threads);
    pop_1 = new_generation_pop_1;
    pop_2 = new_generation_pop_2;

    if (pop_1.empty() || pop_2.empty()) {
      throw "pop magically became empty?";
    }

    if (use_selection) {
      fitness_pop_1 = std::vector<double>(pop_1.size(), 0.0);
      fitness_pop_2 = std::vector<double>(pop_2.size(), 0.0);
      for (size_t i = 0; i < pop_1.size(); ++i) {
        fitness_pop_1[i] = calculate_fitness(pop_1[i], select,
                                             marker_positions,
                                             multiplicative_selection);
      }
      for (size_t i = 0; i < pop_2.size(); ++i) {
        fitness_pop_2[i] = calculate_fitness(pop_2[i], select,
                                             marker_positions,
                                             multiplicative_selection);
      }
      max_fitness_pop_1 = *std::max_element(fitness_pop_1.begin(),
                                            fitness_pop_1.end());
      max_fitness_pop_2 = *std::max_element(fitness_pop_2.begin(),
                                            fitness_pop_2.end());
    }

    if (t % updateFreq == 0 && verbose) {
      Rcout << "**";
    }

    if (t > 2 && is_fixed(pop_1) && is_fixed(pop_2)) {
      if (verbose) Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
      force_output();
      std::vector< std::vector < Fish_emp > > output;
      output.push_back(pop_1);
      output.push_back(pop_2);
      return(output);
    }

    Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcout << "\n";
  std::vector< std::vector< Fish_emp > > output;
  output.push_back(pop_1);
  output.push_back(pop_2);
  return(output);
}

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
  rnd_t rndgen;

  std::vector< Fish_emp > Pop_1;
  std::vector< Fish_emp > Pop_2;
  std::vector<int> founder_labels = {0, 1, 2, 3, 4};
  int number_of_alleles = founder_labels.size();

  std::vector<double> marker_positions(marker_positions_R.begin(),
                                       marker_positions_R.end());

  std::vector<double> track_markers(track_markers_R.begin(),
                                    track_markers_R.end());

  std::vector<size_t> pop_size(2);
  pop_size[0] = static_cast<size_t>(pop_sizes[0]);
  pop_size[1] = static_cast<size_t>(pop_sizes[1]);

  int number_of_markers = track_markers.size();

  emp_genome emp_gen(marker_positions);
  if (static_cast<size_t>(recombination_map.size()) ==
      static_cast<size_t>(marker_positions.size())) {
    std::vector<double> recom_map(recombination_map.begin(),
                                  recombination_map.end());

    emp_gen = emp_genome(recom_map);
    morgan = std::accumulate(recom_map.begin(),
                             recom_map.end(),
                             0.0);
  }

  std::vector< std::vector< double >> substitution_matrix(4,
                                                          std::vector<double>(4));
  if (mutation_rate > 0) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        substitution_matrix[i][j] = substitution_matrix_R(i, j);
      }
    }
  }

  if (input_population_1[0] > -1e4) {

    if (pop_size.size() != 2) {
      stop("pop_size.size() != 2, need two separate population sizes as input");
    }

    Pop_1 = convert_numeric_matrix_to_fish_vector(input_population_1);
    Pop_2 = convert_numeric_matrix_to_fish_vector(input_population_2);

    if (static_cast<size_t>(Pop_1.size()) !=
        static_cast<size_t>(pop_size[0])) {
      //   the populations have to be populated from the parents!
      std::vector< Fish_emp > Pop_1_new(pop_size[0]);
      for(size_t j = 0; j < pop_size[0]; ++j) {
        int index = rndgen.random_number(Pop_1.size());
        Pop_1_new[j] = Pop_1[index];
      }
    }
    if (static_cast<size_t>(Pop_2.size()) !=
        static_cast<size_t>(pop_size[1])) {
      std::vector< Fish_emp > Pop_2_new(pop_size[1]);
      for (int j = 0; j < pop_size[1]; ++j) {
        int index = rndgen.random_number(Pop_2.size());
        Pop_2_new[j] = Pop_2[index];
      }
      Pop_2 = Pop_2_new;
    }
  }

  // 5 columns: time, loc, anc, type, population
  arma::mat frequencies_table(number_of_markers * number_of_alleles * total_runtime * 2, 5);
  arma::mat initial_frequencies = update_all_frequencies_tibble_dual_pop(Pop_1,
                                                                         Pop_2,
                                                                         track_markers,
                                                                         marker_positions,
                                                                         0,
                                                                         morgan);

  std::vector< std::vector< Fish_emp > > output_populations;
  if (verbose) {
    Rcout << "starting simulation\n"; force_output();
  }
  output_populations = simulate_two_populations(Pop_1,
                                                Pop_2,
                                                marker_positions,
                                                select,
                                                pop_size,
                                                total_runtime,
                                                morgan,
                                                verbose,
                                                frequencies_table,
                                                track_frequency,
                                                track_markers,
                                                multiplicative_selection,
                                                number_of_alleles,
                                                founder_labels,
                                                migration_rate,
                                                mutation_rate,
                                                substitution_matrix,
                                                rndgen,
                                                emp_gen,
                                                num_threads);

  if (verbose) { Rcout << "done simulating\n"; force_output(); }

  arma::mat final_frequencies =
    update_all_frequencies_tibble_dual_pop(output_populations[0],
                                           output_populations[1],
                                                             track_markers,
                                                             marker_positions,
                                                             total_runtime,
                                                             morgan);

  std::vector<double> junctions;
  return List::create( Named("population_1") = convert_to_list(output_populations[0],
                             marker_positions),
                             Named("population_2") = convert_to_list(output_populations[1],
                                   marker_positions),
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
