//
//  selection.cpp
//
//
//  Created by Thijs Janzen on 28/02/2018.
//
//
#include <vector>
#include <cstdlib>
#include <numeric>
#include <cmath>

#include <vector>
#include <algorithm>

#include <thread>

#include "Fish.h"
#include "random_functions.h"
#include "helper_functions.h"
#include "util.h"

#include <RcppParallel.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

void update_pop(const std::vector<Fish>& Pop,
                std::vector<Fish>& new_generation,
                size_t pop_size,
                double morgan,
                const std::vector<double>& fitness,
                double maxFitness,
                bool use_selection,
                int num_threads,
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
      while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness, rndgen);
    } else {
      index1 = rndgen.random_number(source_size);
      index2 = rndgen.random_number(source_size);
      while(index2 == index1) index2 = rndgen.random_number(source_size);
    }

    new_generation[i] = mate(Pop[index1], Pop[index2], morgan, rndgen);
  });
}

// Returns the final population; 'generations_run' is set to the number of
// generations simulated (smaller than total_runtime if the population fixed).
std::vector< Fish > simulate_Population(const std::vector< Fish>& sourcePop,
                                        const std::vector< std::array<double, 5>>& select,
                                        size_t pop_size,
                                        int total_runtime,
                                        double morgan,
                                        bool verbose,
                                        arma::mat& frequencies,
                                        bool track_frequency,
                                        const NumericVector& track_markers,
                                        bool track_junctions,
                                        std::vector<double>& junctions,
                                        bool multiplicative_selection,
                                        const std::vector<int>& founder_labels,
                                        rnd_t& rndgen,
                                        int num_threads,
                                        int& generations_run) {

  const bool use_selection = !select.empty();

  std::vector<Fish> Pop = sourcePop;
  std::vector<double> fitness;
  double maxFitness = -1;

  auto update_fitness = [&]() {
    fitness.resize(Pop.size());
    for (size_t i = 0; i < Pop.size(); ++i) {
      fitness[i] = calculate_fitness(Pop[i], select, multiplicative_selection);
    }
    maxFitness = *std::max_element(fitness.begin(), fitness.end());
    if (maxFitness > 0.0 &&
        std::count_if(fitness.begin(), fitness.end(),
                      [](double f) { return f > 0.0; }) < 2) {
      stop("fewer than two individuals with positive fitness remain");
    }
  };

  if (use_selection) update_fitness();

  int updateFreq = total_runtime / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(verbose) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }

  generations_run = 0;
  for(int t = 0; t < total_runtime; ++t) {
    if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));

    if(track_frequency) {
      const int num_labels = founder_labels.size();
      const int time_block = track_markers.size() * num_labels; // number of markers times number of alleles
      for(int i = 0; i < track_markers.size(); ++i) {
        if(track_markers[i] < 0) break;
        arma::mat local_mat = update_frequency_tibble(Pop,
                                                      track_markers[i],
                                                      founder_labels,
                                                      t,
                                                      morgan);

        int start_add_marker = t * time_block + i * num_labels;
        for(int j = 0; j < num_labels; ++j) {
          for(size_t k = 0; k < 4; ++k) {
            frequencies(start_add_marker + j, k)  = local_mat(j, k);
          }
        }
      }
    }

    std::vector<Fish> new_generation(pop_size);

    update_pop(Pop, new_generation, pop_size,
               morgan, fitness, maxFitness, use_selection,
               num_threads, rndgen.next_seed());

    Pop.swap(new_generation);
    generations_run = t + 1;

    if (use_selection) update_fitness();

    if (t % updateFreq == 0 && verbose) {
      Rcout << "**";
    }

    if (t > 2 && is_fixed(Pop)) {
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
    rnd_t rndgen(seed_from_r());

    std::vector< Fish > Pop;
    int number_of_alleles = number_of_founders;
    std::vector<int> founder_labels;

    track_markers = scale_markers(track_markers, morgan);

    if (input_population[0] > -1e4) {
      if (verbose) { Rcout << "found input population, converting\n"; force_output(); }

      Pop = convert_NumericVector_to_fishVector(input_population);

      number_of_founders = 0;
      for (auto it = Pop.begin(); it != Pop.end(); ++it) {
        update_founder_labels((*it).chromosome1, founder_labels);
        update_founder_labels((*it).chromosome2, founder_labels);
      }

      number_of_alleles = founder_labels.size();

      if (Pop.size() != pop_size) {
        // the new population has to be seeded from the input!
        std::vector< Fish > Pop_new;
        for (size_t j = 0; j < pop_size; ++j) {
          int index = rndgen.random_number(Pop.size());
          Pop_new.push_back(Pop[index]);
        }
        Pop = Pop_new;
      }
    } else {
      if (verbose) { Rcout << "starting generating random starting population\n"; force_output(); }
      for (size_t i = 0; i < pop_size; ++i) {
        int founder_1 = draw_random_founder(starting_proportions, rndgen);
        int founder_2 = draw_random_founder(starting_proportions, rndgen);

        Fish p1 = Fish( founder_1 );
        Fish p2 = Fish( founder_2 );

        Pop.push_back(mate(p1,p2, morgan, rndgen));
      }
      for (int i = 0; i < number_of_alleles; ++i) {
        founder_labels.push_back(i);
      }
    }

    arma::mat frequencies_table;

    if (track_frequency) {
      int number_of_markers = track_markers.size();
      arma::mat x(number_of_markers * number_of_alleles * total_runtime, 4); // 4 columns: time, loc, anc, type
      frequencies_table = x;
    }

    arma::mat initial_frequencies = update_all_frequencies_tibble(Pop,
                                                                  track_markers,
                                                                  founder_labels,
                                                                  0,
                                                                  morgan);

    std::vector< std::array<double, 5> > select_matrix =
      convert_select_from_r(select, morgan);

    std::vector<double> junctions;
    int generations_run = 0;
    std::vector<Fish> outputPop = simulate_Population(Pop,
                                                      select_matrix,
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
                                                      founder_labels,
                                                      rndgen,
                                                      num_threads,
                                                      generations_run);

    if (track_frequency) {
      // drop the rows of generations that were not simulated (early fixation)
      const size_t rows_per_generation = track_markers.size() * number_of_alleles;
      frequencies_table.resize(rows_per_generation * generations_run, 4);
    }

    if (verbose) {
      Rcout << "done simulating\n";
    }
    arma::mat final_frequencies = update_all_frequencies_tibble(outputPop,
                                                                track_markers,
                                                                founder_labels,
                                                                generations_run,
                                                                morgan);

    return List::create( Named("population") = convert_to_list(outputPop),
                         Named("frequencies") = frequencies_table,
                         Named("initial_frequencies") = initial_frequencies,
                         Named("final_frequencies") = final_frequencies,
                         Named("junctions") = junctions);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
  return NA_REAL;
}