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

#include "Fish.h"
#include "random_functions.h"
#include "helper_functions.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

template <typename FISH>
std::vector< FISH > simulate_Population(const std::vector<FISH>& sourcePop,
                                        const NumericMatrix& select,
                                        int pop_size,
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
                                        const std::vector<int>& founder_labels) {

  bool use_selection = false;
  if(select(1, 1) >= 0) use_selection = true;

  std::vector<FISH> Pop = sourcePop;
  std::vector<double> fitness;
  double maxFitness = -1;

  if(use_selection) {
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
      double fit = calc_fitness((*it), select, multiplicative_selection);
      if(fit > maxFitness) maxFitness = fit;

      fitness.push_back(fit);
    }
  }

  int updateFreq = total_runtime / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(verbose) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }

  for(int t = 0; t < total_runtime; ++t) {
    if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));

    if(track_frequency) {
      for(int i = 0; i < track_markers.size(); ++i) {
        if(track_markers[i] < 0) break;
        arma::mat local_mat = update_frequency_tibble<FISH>(Pop,
                                                            track_markers[i],
                                                            founder_labels,
                                                            t,
                                                            morgan);

        // now we have to find where to copy local_mat into frequencies
        int time_block = track_markers.size() * founder_labels.size(); // number of markers times number of alleles

        int start_add_time = t * time_block;
        int start_add_marker = i * founder_labels.size() + start_add_time;

        for(size_t j = 0; j < founder_labels.size(); ++j) {
          for(size_t k = 0; k < 4; ++k) {
            frequencies(start_add_marker + j, k)  = local_mat(j, k);
          }
        }
      }
    }

    std::vector<FISH> newGeneration(pop_size);
    std::vector<double> newFitness;
    double newMaxFitness = -1.0;
    for (int i = 0; i < pop_size; ++i)  {
      int index1 = 0;
      int index2 = 0;
      if (use_selection) {
        index1 = draw_prop_fitness(fitness, maxFitness);
        index2 = draw_prop_fitness(fitness, maxFitness);
        while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness);
      } else {
        index1 = random_number( (int)Pop.size() );
        index2 = random_number( (int)Pop.size() );
        while(index2 == index1) index2 = random_number( (int)Pop.size() );
      }

      newGeneration[i] = FISH(Pop[index1].gamete(morgan),
                              Pop[index2].gamete(morgan));
      double fit = -2.0;
      if(use_selection) fit = calc_fitness(newGeneration[i],
         select, multiplicative_selection);
      if(fit > newMaxFitness) newMaxFitness = fit;

      newFitness.push_back(fit);
    }

    if (t % updateFreq == 0 && verbose) {
      Rcout << "**";
    }

    if (t > 2 && is_fixed(Pop) && verbose) {
      Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
      R_FlushConsole();
      return(Pop);
    }

    Rcpp::checkUserInterrupt();

    Pop.swap(newGeneration);
    fitness.swap(newFitness);
    maxFitness = newMaxFitness;
  }
  if(verbose) Rcout << "\n";
  return(Pop);
}

// [[Rcpp::export]]
Rcpp::List simulate_cpp(Rcpp::NumericVector input_population,
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
                        int seed) {

  set_seed(seed);
  set_poisson(morgan);

  typedef Fish<chromosome_junctions> FOCAL_ANIMAL;

  std::vector< FOCAL_ANIMAL > Pop;
  int number_of_alleles = number_of_founders;
  std::vector<int> founder_labels;

  track_markers = scale_markers(track_markers, morgan);

  if (input_population[0] > -1e4) {

    Pop = convert_NumericVector_to_fishVector(input_population);

    number_of_founders = 0;
    for (auto it = Pop.begin(); it != Pop.end(); ++it) {
      update_founder_labels((*it).chromosome1, founder_labels);
      update_founder_labels((*it).chromosome2, founder_labels);
    }
    number_of_alleles = founder_labels.size();

    if (Pop.size() != pop_size) {
      // the new population has to be seeded from the input!
      std::vector< FOCAL_ANIMAL > Pop_new;
      for (size_t j = 0; j < pop_size; ++j) {
        int index = random_number(Pop.size());
        Pop_new.push_back(Pop[index]);
      }
      Pop = Pop_new;
    }
  } else {
    for (size_t i = 0; i < pop_size; ++i) {
      int founder_1 = draw_random_founder(starting_proportions);
      int founder_2 = draw_random_founder(starting_proportions);

      FOCAL_ANIMAL p1 = FOCAL_ANIMAL( founder_1 );
      FOCAL_ANIMAL p2 = FOCAL_ANIMAL( founder_2 );

      Pop.emplace_back( FOCAL_ANIMAL(p1.gamete(morgan),
                                     p2.gamete(morgan)) );
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

  arma::mat initial_frequencies = update_all_frequencies_tibble<FOCAL_ANIMAL>(Pop,
                                                                              track_markers,
                                                                              founder_labels,
                                                                              0,
                                                                              morgan);

  std::vector<double> junctions;
  std::vector<FOCAL_ANIMAL> outputPop =
    simulate_Population< FOCAL_ANIMAL >(Pop,
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
                                        founder_labels);

  arma::mat final_frequencies = update_all_frequencies_tibble<FOCAL_ANIMAL>(outputPop,
                                                                            track_markers,
                                                                            founder_labels,
                                                                            total_runtime,
                                                                            morgan);

  return List::create( Named("population") = convert_to_list(outputPop),
                       Named("frequencies") = frequencies_table,
                       Named("initial_frequencies") = initial_frequencies,
                       Named("final_frequencies") = final_frequencies,
                       Named("junctions") = junctions);
}