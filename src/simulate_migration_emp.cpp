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
                     int &index) {

  Fish_emp parent;
  index = -1;

  if (uniform() < migration_rate) {
    // migration
    if(use_selection) {
      index = draw_prop_fitness(fitness_migr, max_fitness_migr);
    } else {
      index = random_number( (int)pop_2.size() );
    }
    parent = pop_2[index];
    index = index + pop_1.size();
    // to ensure different indices for pop_1 and pop_2
  } else {
    if(use_selection)  {
      index = draw_prop_fitness(fitness_source, max_fitness_source);
    } else {
      index = random_number( (int)pop_1.size() );
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
                                      std::vector< double >& new_fitness,
                                      double size_in_morgan) {

  std::vector<Fish_emp> new_generation(pop_size);

  for (size_t i = 0; i < pop_size; ++i)  {
    int index1, index2;
    Fish_emp parent1 = draw_parent(pop_1, pop_2, migration_rate,
                                   use_selection,
                                   fitness_source, fitness_migr,
                                   max_fitness_source, max_fitness_migr,
                                   index1);
    Fish_emp parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                   use_selection,
                                   fitness_source, fitness_migr,
                                   max_fitness_source, max_fitness_migr,
                                   index2);
    while (index1 == index2) {
      parent2 = draw_parent(pop_1, pop_2, migration_rate,
                            use_selection,
                            fitness_source, fitness_migr,
                            max_fitness_source, max_fitness_migr,
                            index2);
    }

    new_generation[i] = Fish_emp(parent1.gamete(recompos()),
                                 parent2.gamete(recompos()));

    double fit = -2.0;
    if (use_selection) {
        fit = calculate_fitness(new_generation[i], select,
                  marker_positions, multiplicative_selection);

      new_fitness[i] = fit;
    }
  }
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
    double migration_rate) {

  bool use_selection = false;
  if (select(1, 1) >= 0) use_selection = true;

  std::vector<Fish_emp> pop_1 = source_pop_1;
  std::vector<Fish_emp> pop_2 = source_pop_2;

  std::vector<double> fitness_pop_1(pop_1.size(), 0.0);
  std::vector<double> fitness_pop_2(pop_2.size(), 0.0);

  if (use_selection) {
    for (int j = 0; j < select.nrow(); ++j) {
      if (select(j, 4) < 0) break; // these entries are only for tracking, not for selection calculations
      double local_max_fitness = 0.0;
      for (int i = 1; i < 4; ++i) {
        if (select(j, i) > local_max_fitness) {
          local_max_fitness = select(j, i);
        }
      }
    }

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
                                                               new_fitness_pop_1,
                                                               morgan);

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
                                                               new_fitness_pop_2,
                                                               morgan);
    pop_1 = new_generation_pop_1;
    pop_2 = new_generation_pop_2;
    fitness_pop_1 = new_fitness_pop_1;
    fitness_pop_2 = new_fitness_pop_2;
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
                            int seed) {
  set_seed(seed);
  set_poisson(morgan);

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




  if (input_population_1[0] > -1e4) {
    //  if (verbose) Rcout << "Found input populations\n";  force_output();

    Pop_1 = convert_numeric_matrix_to_fish_vector(input_population_1);
    Pop_2 = convert_numeric_matrix_to_fish_vector(input_population_2);

    //  Rcout << "done converting\n"; force_output();

    if (Pop_1.size() != pop_size[0]) {
      //   Rcout << "drawing pop 1: " << pop_size[0] << " from: " << Pop_1.size() << "\n"; force_output();
      // the populations have to be populated from the parents!
      std::vector< Fish_emp > Pop_1_new;
      for(int j = 0; j < pop_size[0]; ++j) {
        int index = random_number(Pop_1.size());
        Pop_1_new.push_back(Pop_1[index]);
      }
      Pop_1 = Pop_1_new;
      Pop_1_new.clear(); // free up memory
    } else {
      Rcpp::stop("can not run without input data");
    }
    //   Rcout << "drawn pop 1\n"; force_output();

    if (Pop_2.size() != pop_size[1]) {
      std::vector< Fish_emp > Pop_2_new;
      //  Rcout << "drawing pop 2: " << pop_size[1] << "\n"; force_output();
      for (int j = 0; j < pop_size[1]; ++j) {
        int index = random_number(Pop_2.size());
        Pop_2_new.push_back(Pop_2[index]);
      }
      Pop_2 = Pop_2_new;
      Pop_2_new.clear(); // free up memory
    }

  }
  // Rcout << "initial frequencies\n"; force_output();
  int number_of_markers = track_markers.size();
  // 5 columns: time, loc, anc, type, population
  arma::mat frequencies_table(number_of_markers * number_of_alleles * total_runtime * 2, 5);
  arma::mat initial_frequencies = update_all_frequencies_tibble_dual_pop(Pop_1,
                                                                         Pop_2,
                                                                         track_markers,
                                                                         marker_positions,
                                                                         0,
                                                                         morgan);

  std::vector<double> junctions;
  std::vector< std::vector< Fish_emp > > output_populations;
  // Rcout << "starting simulation\n"; force_output();
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
                                                migration_rate);

  // Rcout << "done simulating\n"; force_output();

  arma::mat final_frequencies =
    update_all_frequencies_tibble_dual_pop(output_populations[0],
                                           output_populations[1],
                                           track_markers,
                                           marker_positions,
                                           total_runtime,
                                           morgan);

  return List::create( Named("population_1") = convert_to_list(output_populations[0],
                             marker_positions),
                       Named("population_2") = convert_to_list(output_populations[1],
                             marker_positions),
                       Named("frequencies") = frequencies_table,
                       Named("initial_frequencies") = initial_frequencies,
                       Named("final_frequencies") = final_frequencies,
                       Named("junctions") = junctions);
}