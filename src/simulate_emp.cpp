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
#include "util.h"


#include <RcppParallel.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;


void update_pop_emp(const std::vector<Fish_emp>& Pop,
                    std::vector<Fish_emp>& new_generation,
                    size_t pop_size,
                    double morgan,
                    const std::vector<double>& fitness,
                    double maxFitness,
                    bool use_selection,
                    int num_threads,
                    const emp_genome& emp_gen,
                    uint64_t seed) {

  if (new_generation.size() != pop_size) {
    stop("new_generation wrong size");
  }
  const int source_size = static_cast<int>(Pop.size());
  if (source_size < 2) {
    stop("need at least two individuals to produce offspring");
  }

  parallel_blocks(pop_size, num_threads, seed, [&](size_t i, rnd_t& rndgen) {
    int index1 = 0;
    int index2 = 0;
    if (use_selection) {
      index1 = draw_prop_fitness(fitness, maxFitness, rndgen);
      index2 = draw_prop_fitness(fitness, maxFitness, rndgen);
      while(index1 == index2) {
        index2 = draw_prop_fitness(fitness, maxFitness, rndgen);
      }
    } else {
      index1 = rndgen.random_number(source_size);
      index2 = rndgen.random_number(source_size);
      while(index1 == index2) {
        index2 = rndgen.random_number(source_size);
      }
    }

    new_generation[i] = Fish_emp(Pop[index1].gamete(morgan, rndgen, emp_gen),
                                 Pop[index2].gamete(morgan, rndgen, emp_gen));
  });
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
                                                const std::vector<double>& track_markers,
                                                bool multiplicative_selection,
                                                double mutation_rate,
                                                const std::vector<std::vector<double>>& sub_matrix,
                                                rnd_t& rndgen,
                                                const emp_genome& emp_gen,
                                                int num_threads,
                                                int& generations_run) {

  size_t num_alleles = 5;
  const bool use_selection = !select_matrix.empty();

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

  generations_run = 0;
  for(int t = 0; t < total_runtime; ++t) {
    // Rcout << t << " " << Pop.size() << "\n"; force_output();
    if(track_frequency) {
      // Rcout << t << " update frequency tibble\n";
      for(size_t i = 0; i < track_markers.size(); ++i) {
        if(track_markers[i] < 0) break;

        int index = find_location(marker_positions, track_markers[i]);
        if (index < 0) {
          throw std::runtime_error("tracked marker not found in data");
        }

        std::vector<std::vector<double>> local_mat = update_frequency_tibble(Pop,
                                                                             index,
                                                                             track_markers[i],
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
                   emp_gen,
                   rndgen.next_seed());

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

    Pop.swap(newGeneration);
    generations_run = t + 1;

    // with mutation, a fixed population can still change
    if (mutation_rate <= 0 && t > 2 && is_fixed(Pop)) {
      if (verbose) Rcout << "\n After " << generations_run << " generations, the population has become completely homozygous and fixed\n";
      R_FlushConsole();
      return(Pop);
    }

    Rcpp::checkUserInterrupt();
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
    rnd_t rndgen(seed_from_r());

    std::vector< std::array<double, 5> > select_matrix =
      convert_select_from_r(select);

    std::vector<double> marker_positions(marker_positions_R.begin(),
                                         marker_positions_R.end());

    std::vector<double> track_markers(track_markers_R.begin(),
                                      track_markers_R.end());

    // default: crossovers proportional to the distance between markers
    emp_genome emp_gen = emp_genome::from_positions(marker_positions);
    if (recombination_map.size() > 0 &&
        static_cast<size_t>(recombination_map.size()) ==
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

    int generations_run = 0;
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
                                                               num_threads,
                                                               generations_run);

    if (track_frequency) {
      // drop the rows of generations that were not simulated (early fixation)
      const size_t rows_per_generation = track_markers.size() * number_of_alleles;
      frequencies_table.resize(rows_per_generation * generations_run, 4);
    }

    arma::mat final_frequencies = update_all_frequencies_tibble(output_pop,
                                                                marker_positions,
                                                                marker_positions,
                                                                generations_run,
                                                                morgan);

    return List::create( Named("population") = convert_to_list(output_pop,
                               marker_positions),
                               Named("frequencies") = frequencies_table,
                               Named("initial_frequencies") = initial_frequencies,
                               Named("final_frequencies") = final_frequencies);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
  return NA_REAL;
}