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
                                      const NumericMatrix& substitution_matrix,
                                      rnd_t& rndgen,
                                      const emp_genome& emp_gen,
                                      int num_threads) {

  std::vector<Fish_emp> new_generation(pop_size);

  int seed_index = 0;
  std::mutex mutex;
  int num_seeds = num_threads * 2; // tbb might re-start threads due to the load-balancer
  if (num_threads == -1) {
    num_seeds = 20;
  }
  std::vector< int > seed_values(num_seeds);

  for (int i = 0; i < num_seeds; ++i) {
    seed_values[i] = rndgen.random_number(INT_MAX); // large value
  }

  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

  tbb::parallel_for(
    tbb::blocked_range<unsigned>(0, pop_size),
    [&](const tbb::blocked_range<unsigned>& r) {

      rnd_t rndgen2(seed_values[seed_index]);
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


        int index1, index2;
        Fish_emp parent1 = draw_parent(pop_1, pop_2, migration_rate,
                                       use_selection,
                                       fitness_source, fitness_migr,
                                       max_fitness_source, max_fitness_migr,
                                       index1,
                                       rndgen);
        Fish_emp parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                       use_selection,
                                       fitness_source, fitness_migr,
                                       max_fitness_source, max_fitness_migr,
                                       index2,
                                       rndgen);
        while (index1 == index2) {
          parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                use_selection,
                                fitness_source, fitness_migr,
                                max_fitness_source, max_fitness_migr,
                                index2,
                                rndgen);
        }

        new_generation[i] = Fish_emp(parent1.gamete(size_in_morgan, rndgen, emp_gen),
                                     parent2.gamete(size_in_morgan, rndgen, emp_gen));

        if (mutation_rate > 0)
          mutate(new_generation[i], substitution_matrix, mutation_rate, rndgen);

        double fit = -2.0;
      }
    });
    return new_generation;
}

std::vector< std::vector< Fish_emp > > simulate_two_populations(
    const std::vector< Fish_emp>& source_pop_1,
    const std::vector< Fish_emp>& source_pop_2,
    const std::vector<double>& marker_positions,
    const NumericMatrix& select,
    const NumericVector& pop_size,
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
    const NumericMatrix& substitution_matrix,
    rnd_t& rndgen,
    const emp_genome& emp_gen,
    int num_threads) {

  bool use_selection = false;
  if (select(1, 1) >= 0) use_selection = true;

  std::vector<Fish_emp> pop_1 = source_pop_1;
  std::vector<Fish_emp> pop_2 = source_pop_2;

  std::vector<double> fitness_pop_1(pop_1.size(), 0.0);
  std::vector<double> fitness_pop_2(pop_2.size(), 0.0);

  if (use_selection) {
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

    std::vector<double> new_fitness_pop_1(pop_1.size(), 0.0);
    std::vector<double> new_fitness_pop_2(pop_2.size(), 0.0);

    std::vector<Fish_emp> new_generation_pop_1 = next_pop_migr(pop_1, // resident
                                                               pop_2, // migrants
                                                               marker_positions,
                                                               pop_size[0],
                                                                       fitness_pop_1,
                                                                       fitness_pop_2,
                                                                       max_fitness_pop_1,
                                                                       max_fitness_pop_2,
                                                                       select,
                                                                       use_selection,
                                                                       multiplicative_selection,
                                                                       migration_rate,
                                                                       morgan,
                                                                       mutation_rate,
                                                                       substitution_matrix,
                                                                       rndgen,
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
                                                                       select,
                                                                       use_selection,
                                                                       multiplicative_selection,
                                                                       migration_rate,
                                                                       morgan,
                                                                       mutation_rate,
                                                                       substitution_matrix,
                                                                       rndgen,
                                                                       emp_gen,
                                                                       num_threads);
    pop_1 = new_generation_pop_1;
    pop_2 = new_generation_pop_2;

    if (use_selection) {
      for (int i = 0; i < pop_1.size(); ++i) {
        fitness_pop_1[i] = calculate_fitness(pop_1[i], select,
                              marker_positions, multiplicative_selection);
      }
      for (int i = 0; i < pop_2.size(); ++i) {
        fitness_pop_2[i] = calculate_fitness(pop_2[i], select,
                                             marker_positions, multiplicative_selection);
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
                                NumericVector pop_size,
                                int total_runtime,
                                double morgan,
                                bool verbose,
                                bool track_frequency,
                                const NumericVector& track_markers_R,
                                bool multiplicative_selection,
                                double migration_rate,
                                int seed,
                                double mutation_rate,
                                const NumericMatrix& substitution_matrix,
                                int num_threads) {

  rnd_t rndgen(seed);

  std::vector< Fish_emp > Pop_1;
  std::vector< Fish_emp > Pop_2;
  std::vector<int> founder_labels = {0, 1, 2, 3, 4};
  int number_of_alleles = founder_labels.size();

  std::vector<double> marker_positions(marker_positions_R.begin(),
                                       marker_positions_R.end());

  auto inv_max_marker_pos = 1.0 / (*std::max_element(marker_positions.begin(),
                                                     marker_positions.end()));

  for (auto& i : marker_positions) {
    i *= inv_max_marker_pos;
  }

  std::vector<double> track_markers(track_markers_R.begin(),
                                    track_markers_R.end());

  int number_of_markers = track_markers.size();
  if (verbose) Rcout << number_of_markers << "\n";

  if (*std::max_element(track_markers_R.begin(), track_markers_R.end()) > 1) {
    for (auto& i : track_markers) {
      i *= inv_max_marker_pos;
      if (verbose) {
        Rcout << i << " "; force_output();
      }
    }
  }

  if (verbose) { Rcout << "markers loaded\n"; force_output(); }

  if (select.nrow() > 0) {
    if (select(0, 0) > 10) {// location in bp
      for (int i = 0; i < select.nrow(); ++i) {
        select(i, 0) = select(i, 0) * inv_max_marker_pos;
      }
    } else {
      if (select(0, 0) > 1) {
        for (int i = 0; i < select.nrow(); ++i) {
          select(i, 0) = select(i, 0) / morgan;
        }
      }
    }
  }


  emp_genome emp_gen(marker_positions);

  if (input_population_1[0] > -1e4) {
    if (verbose) { Rcout << "Found input populations\n";  force_output(); }

    if (pop_size.size() != 2) {
      stop("pop_size.size() != 2, need two separate population sizes as input");
    }

    Pop_1 = convert_numeric_matrix_to_fish_vector(input_population_1);
    Pop_2 = convert_numeric_matrix_to_fish_vector(input_population_2);

    if (verbose) { Rcout << "done converting\n"; force_output(); }
    if (verbose) { Rcout << "pop1: " << Pop_1.size() << "\n"; force_output(); }
    if (verbose) { Rcout << "pop2: " << Pop_2.size() << "\n"; force_output(); }
    if (verbose) { Rcout << pop_size[0] << " " << pop_size[1] << "\n"; force_output(); }

    if (static_cast<size_t>(Pop_1.size()) !=
        static_cast<size_t>(pop_size[0])) {
      if (verbose) {Rcout << "drawing pop 1: " << pop_size[0] <<
        " from: " << Pop_1.size() << "\n"; force_output(); }
    //   the populations have to be populated from the parents!
      std::vector< Fish_emp > Pop_1_new(pop_size[0]);
      for(size_t j = 0; j < pop_size[0]; ++j) {
        int index = random_number(Pop_1.size());
        Pop_1_new[j] = Pop_1[index];
    //  Rcout << "done converting\n"; force_output();

    if (verbose)  {Rcout << "drawn pop 1\n"; force_output();}

    if (static_cast<size_t>(Pop_2.size()) !=
        static_cast<size_t>(pop_size[1])) {
      std::vector< Fish_emp > Pop_2_new(pop_size[1]);
      if (verbose) { Rcout << "drawing pop 2: " << pop_size[1] << " " <<
                              Pop_2.size() << "\n"; force_output(); }
      for (int j = 0; j < pop_size[1]; ++j) {
        int index = random_number(Pop_2.size());
        Pop_2_new[j] = Pop_2[index];
      }
      Pop_2 = Pop_2_new;
    }

    if (verbose) {Rcout << "drawing done: " << Pop_1.size() << " " <<
                                              Pop_2.size() << "\n"; force_output();}
  }

  if (verbose) {Rcout << "input_data loaded\n"; force_output();}

  if (verbose) {Rcout << "initial frequencies\n"; force_output();}

  if (verbose) {Rcout << track_markers.size() << "\n"; force_output; }
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
}
