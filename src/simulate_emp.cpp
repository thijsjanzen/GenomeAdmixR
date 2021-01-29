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
#include <tuple>

#include "Fish_emp.h"
#include "random_functions.h"
#include "helper_functions.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;


std::vector< Fish_emp > simulate_population_emp(const std::vector< Fish_emp>& sourcePop,
                                                const NumericMatrix& select_matrix,
                                                const std::vector<double>& marker_positions,
                                                int pop_size,
                                                int total_runtime,
                                                double morgan,
                                                bool verbose,
                                                arma::mat& frequencies,
                                                bool track_frequency,
                                                const std::vector<double>& track_markers,
                                                bool multiplicative_selection) {

  int num_alleles = 5;
  bool use_selection = false;
  if(select_matrix(0, 5) >= 0) use_selection = true;
 /* for(int i = 0; i < select_matrix.nrow(); ++i) {
    for(int j = 0; j < select_matrix.ncol(); ++j) {
      Rcout << select_matrix(i, j) << " ";
    }
    Rcout << "\n";
 }*/


  std::vector<Fish_emp> Pop = sourcePop;
  std::vector<double> fitness;


  if(use_selection) {
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
      double fit = calculate_fitness((*it),
                                     select_matrix,
                                     marker_positions,
                                     multiplicative_selection);

      fitness.push_back(fit);
    }
  }
  double maxFitness = *std::max_element(fitness.begin(), fitness.end());

  int updateFreq = total_runtime / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(verbose) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }

  for(int t = 0; t < total_runtime; ++t) {
  //  Rcout << t << " " << Pop.size() << "\n";
    if(track_frequency) {
 //     Rcout << t << " update frequency tibble\n";
      for(int i = 0; i < track_markers.size(); ++i) {
        if(track_markers[i] < 0) break;

        int index = find_location(marker_positions, track_markers[i]);
        if (index < 0)
          continue;

        std::vector<std::vector<double>> local_mat = update_frequency_tibble(Pop,
                                                                             index,
                                                                             track_markers[i] * morgan,
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
      //Rcout << "tibble updated\n";
    }

    std::vector<Fish_emp> newGeneration(pop_size);
    std::vector<double> newFitness;

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

      newGeneration[i] = Fish_emp(Pop[index1].gamete(recompos()),
                                  Pop[index2].gamete(recompos()));



      if(use_selection) {
          double fit = calculate_fitness(newGeneration[i],
                                                select_matrix,
                                                marker_positions,
                                                multiplicative_selection);

          newFitness.push_back(fit);
      }
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
    maxFitness = *std::max_element(fitness.begin(), fitness.end());
  }
  if(verbose) Rcout << "\n";
  return(Pop);
}

// [[Rcpp::export]]
List simulate_emp_cpp(Rcpp::NumericMatrix input_population,
                      Rcpp::NumericVector marker_positions_R,
                      Rcpp::NumericMatrix select,
                      size_t pop_size,
                      size_t total_runtime,
                      double morgan,
                      bool verbose,
                      bool track_frequency,
                      Rcpp::NumericVector track_markers_R,
                      bool multiplicative_selection,
                      int seed) {

  set_seed(seed);
  set_poisson(morgan);

  std::vector<double> marker_positions(marker_positions_R.begin(),
                                       marker_positions_R.end());

  auto inv_max_marker_pos = 1.0 / (*std::max_element(marker_positions.begin(),
                                                    marker_positions.end()));

//  Rcout << inv_max_marker_pos << "\n"; force_output();

  for (auto& i : marker_positions) {
    i *= inv_max_marker_pos;
  }

  std::vector<double> track_markers(track_markers_R.begin(),
                                    track_markers_R.end());

  if (*std::max_element(track_markers_R.begin(), track_markers_R.end()) > 1) {
    for (auto& i : track_markers) {
      i *= inv_max_marker_pos;
    }
  }

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


  fill_cum_marker_dist(marker_positions);

  std::vector< Fish_emp > Pop;

  std::vector<int> founder_labels = {0, 1, 2, 3, 4};
  int number_of_alleles = founder_labels.size();

  track_markers = scale_markers(track_markers, morgan);

  if (input_population[0] > -1e4) {
 //   Rcout << "converting\n"; force_output();
    Pop = convert_numeric_matrix_to_fish_vector(input_population);

    // the new population has to be seeded from the input!
    std::vector< Fish_emp > Pop_new;
    for (size_t j = 0; j < pop_size; ++j) {
      int index = random_number(Pop.size());
      Pop_new.push_back(Pop[index]);
    }
    Pop = Pop_new;
  } else {
    Rcpp::stop("can not run without input data");
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

 // Rcout << "simulate\n"; force_output();
  std::vector<Fish_emp> output_pop = simulate_population_emp(Pop,
                                                             select,
                                                             marker_positions,
                                                             pop_size,
                                                             total_runtime,
                                                             morgan,
                                                             verbose,
                                                             frequencies_table,
                                                             track_frequency,
                                                             track_markers,
                                                             multiplicative_selection);

 // Rcout << "final frequencies\n"; force_output();
  arma::mat final_frequencies = update_all_frequencies_tibble(output_pop,
                                                              marker_positions,
                                                              marker_positions,
                                                              total_runtime,
                                                              morgan);

  //Rcout << "convert to list\n"; force_output();
  return List::create( Named("population") = convert_to_list(output_pop,
                             marker_positions),
                             Named("frequencies") = frequencies_table,
                             Named("initial_frequencies") = initial_frequencies,
                             Named("final_frequencies") = final_frequencies);
}