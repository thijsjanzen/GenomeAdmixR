//
//  selection.cpp
//
//
//  Created by Thijs Janzen on 28/02/2018.
//
//
#include <iostream>
#include <vector>
#include <fstream>
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

#include <tbb/tbb.h>

void     update_pop(const std::vector<Fish>& Pop,
                    const std::vector<double> fitness,
                    std::vector<Fish>& newGeneration,
                    std::vector<double>& newFitness,
                    const NumericMatrix& select,
                    double maxFitness,
                    double morgan,
                    bool multiplicative_selection,
                    bool use_selection,
                    int num_threads) {


  if (num_threads == 1) {
    for (int i = 0; i < Pop.size(); ++i)  {
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

      newGeneration[i] = mate(Pop[index1], Pop[index2], morgan);

      double fit = -2.0;
      if(use_selection) fit = calculate_fitness(newGeneration[i], select, multiplicative_selection);

      newFitness[i] = fit;
    }
  }

  if (num_threads > 1) {

    int pop_size = Pop.size();

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, pop_size),
       [&](const tbb::blocked_range<unsigned>& r) {

         std::random_device rd;
         std::mt19937 local_rndgen_(rd());

         // initialize randomizers.
         std::uniform_int_distribution<> pop_draw = std::uniform_int_distribution<> (0, Pop.size() - 1);

         for (unsigned i = r.begin(); i < r.end(); ++i) {
           int index1 = 0;
           int index2 = 0;
           if (use_selection) {
             index1 = draw_prop_fitness(fitness, maxFitness, local_rndgen_);
             index2 = draw_prop_fitness(fitness, maxFitness, local_rndgen_);
             while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness, local_rndgen_);
           } else {
             index1 = pop_draw(local_rndgen_);
             index2 = pop_draw(local_rndgen_);
             while(index2 == index1) index2 = pop_draw(local_rndgen_);
           }

           newGeneration[i] = mate_threaded(Pop[index1], Pop[index2], morgan,
                                            local_rndgen_);

           double fit = -2.0;
           if(use_selection) fit = calculate_fitness(newGeneration[i], select, multiplicative_selection);

           newFitness[i] = fit;
         }
       }
      );
  }
  return;
}


std::vector< Fish > simulate_Population(const std::vector< Fish>& sourcePop,
                                        const NumericMatrix& select,
                                        int pop_size,
                                        int total_runtime,
                                        double morgan,
                                        bool progress_bar,
                                        arma::mat& frequencies,
                                        bool track_frequency,
                                        const NumericVector& track_markers,
                                        bool track_junctions,
                                        std::vector<double>& junctions,
                                        bool multiplicative_selection,
                                        int num_alleles,
                                        const std::vector<int>& founder_labels,
                                        int num_threads) {

  //Rcout << "simulate_population: " << multiplicative_selection << "\n";

  bool use_selection = false;
  if(select(1, 1) >= 0) use_selection = true;

  std::vector<Fish> Pop = sourcePop;
  std::vector<double> fitness;
  double maxFitness = -1;

  if(use_selection) {
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
      double fit = calculate_fitness((*it), select, multiplicative_selection);
      if(fit > maxFitness) maxFitness = fit;

      fitness.push_back(fit);
    }
  }

  int updateFreq = total_runtime / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(progress_bar) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }

  tbb::task_scheduler_init _tbb((num_threads > 1) ? num_threads : tbb::task_scheduler_init::automatic);

  for(int t = 0; t < total_runtime; ++t) {

    if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));

    if(track_frequency) {
      for(int i = 0; i < track_markers.size(); ++i) {
        if(track_markers[i] < 0) break;
        arma::mat local_mat = update_frequency_tibble(Pop,
                                                      track_markers[i],
                                                                   founder_labels,
                                                                   t,
                                                                   morgan);

        // now we have to find where to copy local_mat into frequencies
        int time_block = track_markers.size() * founder_labels.size(); // number of markers times number of alleles

        int start_add_time = t * time_block;
        int start_add_marker = i * founder_labels.size() + start_add_time;

        for(int j = 0; j < founder_labels.size(); ++j) {
          for(int k = 0; k < 4; ++k) {
            frequencies(start_add_marker + j, k)  = local_mat(j, k);
          }
        }
      }
    }

    std::vector<Fish> newGeneration(pop_size);
    std::vector<double> newFitness(pop_size);
    double newMaxFitness = -1.0;

    update_pop(Pop,
               fitness,
               newGeneration,
               newFitness,
               select,
               maxFitness,
               morgan,
               multiplicative_selection,
               use_selection,
               num_threads);

    if (t % updateFreq == 0 && progress_bar) {
      Rcout << "**";
    }

    if (t > 2 && is_fixed(Pop)) {
      Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
      R_FlushConsole();
      return(Pop);
    }

    Rcpp::checkUserInterrupt();

    Pop.swap(newGeneration);
    fitness.swap(newFitness);
    maxFitness = *std::max_element(fitness.begin(), fitness.end());
  }
  if(progress_bar) Rcout << "\n";
  return(Pop);
}

// [[Rcpp::export]]
List simulate_cpp(Rcpp::NumericVector input_population,
                  NumericMatrix select,
                  int pop_size,
                  int number_of_founders,
                  Rcpp::NumericVector starting_proportions,
                  int total_runtime,
                  double morgan,
                  bool progress_bar,
                  bool track_frequency,
                  NumericVector track_markers,
                  bool track_junctions,
                  bool multiplicative_selection,
                  int seed,
                  int num_threads) {

  set_seed(seed);
  set_poisson(morgan);

  std::vector< Fish > Pop;
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
      std::vector< Fish > Pop_new;
      for (int j = 0; j < pop_size; ++j) {
        int index = random_number(Pop.size());
        Pop_new.push_back(Pop[index]);
      }
      std::swap(Pop, Pop_new);
    }
  } else {
    for (int i = 0; i < pop_size; ++i) {
      int founder_1 = draw_random_founder(starting_proportions);
      int founder_2 = draw_random_founder(starting_proportions);

      Fish p1 = Fish( founder_1 );
      Fish p2 = Fish( founder_2 );

      Pop.push_back(mate(p1,p2, morgan));
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

  std::vector<double> junctions;
  std::vector<Fish> outputPop = simulate_Population(Pop,
                                                    select,
                                                    pop_size,
                                                    total_runtime,
                                                    morgan,
                                                    progress_bar,
                                                    frequencies_table,
                                                    track_frequency,
                                                    track_markers,
                                                    track_junctions,
                                                    junctions,
                                                    multiplicative_selection,
                                                    number_of_alleles,
                                                    founder_labels,
                                                    num_threads);

  arma::mat final_frequencies = update_all_frequencies_tibble(outputPop,
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