//
//  selection.cpp
//
//
//  Created by Thijs Janzen on 28/02/2018.
//
//
#include <array>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <cmath>

#include <vector>
#include <algorithm>
#include <tuple>
#include <mutex>


#include "Fish_emp.h"
#include "random_functions.h"
#include "helper_functions.h"


#include <RcppParallel.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;


void update_pop_emp(const std::vector<Fish_emp>& Pop,
                    std::vector<Fish_emp>& new_generation,
                    size_t pop_size,
                    double morgan,
                    const std::vector<double>& fitness,
                    const double& maxFitness,
                    bool use_selection,
                    int num_threads,
                    const emp_genome& emp_gen_input) {

  if (Pop.size() != pop_size) {
    stop("wrong size pop");
  }
  if (new_generation.size() != pop_size) {
    stop("new_generation wrong size");
  }

  if (num_threads == 1) {
    rnd_t rndgen2;
    emp_genome local_emp_genome = emp_gen_input;
    for (size_t i = 0; i < pop_size; ++i) {
      int index1 = 0;
      int index2 = 0;
      if (use_selection) {
        index1 = draw_prop_fitness(fitness, maxFitness, rndgen2);
        index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
        while(index1 == index2) {
          index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
        }
      } else {
        index1 = rndgen2.random_number( pop_size );
        index2 = rndgen2.random_number( pop_size );
        while(index1 == index2) {
          index2 = rndgen2.random_number( pop_size );
        }
      }

      new_generation[i] = Fish_emp(Pop[index1].gamete(morgan,
                                                      rndgen2,
                                                      local_emp_genome),
                                                      Pop[index2].gamete(morgan,
                                                                         rndgen2,
                                                                         local_emp_genome));
    }
  } else {

    int seed_index = 0;
    std::mutex mutex;
    int num_seeds = num_threads * 2; // tbb might re-start threads due to the load-balancer
    if (num_threads == -1) {
      num_seeds = 20;
    }

    std::vector< int > seed_values(num_seeds);

    { // ensure that rndgen only does this, and doesn't live outside scope.
      rnd_t rndgen;
      for (int i = 0; i < num_seeds; ++i) {
        seed_values[i] = rndgen.random_number(INT_MAX); // large value
      }
    }

    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, pop_size),
      [&](const tbb::blocked_range<unsigned>& r) {

        emp_genome local_emp_genome;
        rnd_t rndgen2(seed_values[seed_index]);
        {
          std::lock_guard<std::mutex> m(mutex);
          local_emp_genome = emp_gen_input;
          seed_index++;
          if (seed_index >= num_seeds) { // just in case.
            for (int i = 0; i < num_seeds; ++i) {
              seed_values[i] = rndgen2.random_number(INT_MAX);
            }
            seed_index = 0;
          }
        }

        for (unsigned i = r.begin(); i < r.end(); ++i) {
          int index1 = 0;
          int index2 = 0;
          if (use_selection) {
            index1 = draw_prop_fitness(fitness, maxFitness, rndgen2);
            index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
            while(index1 == index2) {
              index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
            }
          } else {
            index1 = rndgen2.random_number( pop_size );
            index2 = rndgen2.random_number( pop_size );
            while(index1 == index2) {
              index2 = rndgen2.random_number( pop_size );
            }
          }

          new_generation[i] = Fish_emp(Pop[index1].gamete(morgan,
                                                          rndgen2,
                                                          local_emp_genome),
                                                          Pop[index2].gamete(morgan,
                                                                             rndgen2,
                                                                             local_emp_genome));
        }
      }
    );
  }

  return;
}



std::vector< Fish_emp > simulate_population_emp(const std::vector< Fish_emp>& sourcePop,
                                                const std::vector<std::array<double, 5>>& select_matrix,
                                                const std::vector<double>& marker_positions,
                                                size_t pop_size,
                                                int total_runtime,
                                                double morgan,
                                                bool verbose,
                                                arma::mat& frequencies,
                                                bool track_frequency,
                                                const std::vector<int>& track_markers,
                                                bool multiplicative_selection,
                                                double mutation_rate,
                                                const std::vector<std::vector<double>>& sub_matrix,
                                                rnd_t& rndgen,
                                                const emp_genome& emp_gen,
                                                int num_threads) {

  size_t num_alleles = 5;
  bool use_selection = false;
  if(select_matrix[0][0] >= 0) use_selection = true;

  std::vector<Fish_emp> Pop = sourcePop;
  std::vector<double> fitness;
  double maxFitness = 0.0;

  if(use_selection) {
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
      double fit = calculate_fitness((*it),
                                     select_matrix,
                                     marker_positions,
                                     multiplicative_selection);

      fitness.push_back(fit);
    }
    maxFitness = *std::max_element(fitness.begin(), fitness.end());
  }

  int updateFreq = total_runtime / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(verbose) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }

  for(int t = 0; t < total_runtime; ++t) {
    // Rcout << t << " " << Pop.size() << "\n"; force_output();
    if(track_frequency) {
      // Rcout << t << " update frequency tibble\n";
      for(size_t i = 0; i < track_markers.size(); ++i) {
        if(track_markers[i] < 0) break;

        int index = find_location(marker_positions, track_markers[i]);
        if (index < 0)
          continue;

        //  Rcout << track_markers[i] << " " << index << "\n"; force_output();
        int pos_in_bp = track_markers[i];
        std::vector<std::vector<double>> local_mat = update_frequency_tibble(Pop,
                                                                             index,
                                                                             pos_in_bp,
                                                                             t);

        // now we have to find where to copy local_mat into frequencies
        int time_block = track_markers.size() * num_alleles; // number of markers times number of alleles

        int start_add_time = t * time_block;
        int start_add_marker = i * num_alleles + start_add_time;

        for(size_t j = 0; j < num_alleles; ++j) {
          for(size_t k = 0; k < 4; ++k) { // 4 columns
            frequencies(start_add_marker + j, k)  = local_mat[j][k];
          }
        }
      }
      //   Rcout << "tibble updated\n"; force_output();
    }

    std::vector<Fish_emp> newGeneration(pop_size);

    update_pop_emp(Pop,
                   newGeneration,
                   pop_size,
                   morgan,
                   fitness,
                   maxFitness,
                   use_selection,
                   num_threads,
                   emp_gen);

    if (mutation_rate > 0) {
      for(auto it = newGeneration.begin(); it != newGeneration.end(); ++it) {
        mutate((*it), sub_matrix, mutation_rate, rndgen);
      }
    }


    if(use_selection) {
      for(size_t i = 0; i < newGeneration.size(); ++i) {
        fitness[i] = calculate_fitness(newGeneration[i],
                                       select_matrix,
                                       marker_positions,
                                       multiplicative_selection);
      }
      maxFitness = *std::max_element(fitness.begin(), fitness.end());
    }

    if (t % updateFreq == 0 && verbose) {
      Rcout << "**";
    }

    if (t > 2 && is_fixed(Pop)) {
      if (verbose) Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
      R_FlushConsole();
      return(Pop);
    }

    Rcpp::checkUserInterrupt();

    Pop.swap(newGeneration);
  }
  if(verbose) Rcout << "\n";
  return(Pop);
}

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
    rnd_t rndgen;

    std::vector< std::array<double, 5> > select_matrix;
    for (size_t i = 0; i < select.nrow(); ++i) {
      std::array<double, 5> row_entry;
      for (size_t j = 0; j < select.ncol(); ++j) {
        row_entry[j] = select(i, j);
      }
      select_matrix.push_back(row_entry);
    }

    std::vector<double> marker_positions(marker_positions_R.begin(),
                                         marker_positions_R.end());

    std::vector<int> track_markers(track_markers_R.begin(),
                                   track_markers_R.end());

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

    std::vector< Fish_emp > Pop;

    std::vector<int> founder_labels = {0, 1, 2, 3, 4};
    int number_of_alleles = founder_labels.size();

    //  track_markers = scale_markers(track_markers, morgan);

    if (input_population(0, 0) > -1e4) {
      Pop = convert_numeric_matrix_to_fish_vector(input_population);

      // the new population has to be seeded from the input!
      std::vector< Fish_emp > Pop_new;
      for (size_t j = 0; j < pop_size; ++j) {
        int index = rndgen.random_number(Pop.size());
        Pop_new.push_back(Pop[index]);
      }
      Pop = Pop_new;
    } else {
      Rcpp::stop("can not run without input data");
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

    arma::mat frequencies_table;

    if (track_frequency) {
      int number_of_markers = track_markers.size();
      arma::mat x(number_of_markers * number_of_alleles * total_runtime, 4); // 4 columns: time, loc, anc, type
      frequencies_table = x;
    }

    arma::mat initial_frequencies = update_all_frequencies_tibble(Pop,
                                                                  marker_positions,
                                                                  marker_positions,
                                                                  0,
                                                                  morgan);

    std::vector<Fish_emp> output_pop = simulate_population_emp(Pop,
                                                               select_matrix,
                                                               marker_positions,
                                                               pop_size,
                                                               total_runtime,
                                                               morgan,
                                                               verbose,
                                                               frequencies_table,
                                                               track_frequency,
                                                               track_markers,
                                                               multiplicative_selection,
                                                               mutation_rate,
                                                               substitution_matrix,
                                                               rndgen,
                                                               emp_gen,
                                                               num_threads);

    arma::mat final_frequencies = update_all_frequencies_tibble(output_pop,
                                                                marker_positions,
                                                                marker_positions,
                                                                total_runtime,
                                                                morgan);

    return List::create( Named("population") = convert_to_list(output_pop,
                               marker_positions),
                               Named("frequencies") = frequencies_table,
                               Named("initial_frequencies") = initial_frequencies,
                               Named("final_frequencies") = final_frequencies);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}