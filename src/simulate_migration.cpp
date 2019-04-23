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
#include <assert.h>

#include <vector>
#include <algorithm>

#include "Fish.h"
#include "random_functions.h"
#include "helper_functions.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

Fish draw_parent(const std::vector< Fish>& pop_1,
                 const std::vector< Fish>& pop_2,
                 double migration_rate,
                 bool use_selection,
                 std::vector< double > fitness_source,
                 std::vector< double > fitness_migr,
                 double max_fitness_source,
                 double max_fitness_migr) {

    Fish parent;
  //  Rcout << "start draw_parent\n";

    if(uniform() < migration_rate) {
        // migration
        int index;

        if(use_selection) {
            index = draw_prop_fitness(fitness_migr, max_fitness_migr);
        } else {
            index = random_number( (int)pop_2.size() );
        }
        assert(index < pop_2.size());
        parent = pop_2[index];
    } else {
        int index;
        if(use_selection)  {
            index = draw_prop_fitness(fitness_source, max_fitness_source);
        } else {
            index = random_number( (int)pop_1.size() );
        }
        assert(index < pop_1.size());
        parent = pop_1[index];
    }
    //Rcout << "end draw_parent\n";
    return(parent);
}



std::vector< Fish > next_pop_migr(const std::vector< Fish>& pop_1,
                                  const std::vector< Fish>& pop_2,
                                  int pop_size,
                                  std::vector< double > fitness_source,
                                  std::vector< double > fitness_migr,
                                  double max_fitness_source,
                                  double max_fitness_migr,
                                  const NumericMatrix& select,
                                  bool use_selection,
                                  bool multiplicative_selection,
                                  double migration_rate,
                                  std::vector< double >& new_fitness,
                                  double& new_max_fitness,
                                  double size_in_morgan) {

    std::vector<Fish> new_generation;
    new_fitness.clear();
    new_max_fitness = -1.0;
    for(int i = 0; i < pop_size; ++i)  {
        Fish parent1 = draw_parent(pop_1, pop_2, migration_rate,
                                   use_selection,
                                   fitness_source, fitness_migr,
                                   max_fitness_source, max_fitness_migr);
        Fish parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                   use_selection,
                                   fitness_source, fitness_migr,
                                   max_fitness_source, max_fitness_migr);

        while(parent1 == parent2) {
            parent2 = draw_parent(pop_1, pop_2, migration_rate,
                                  use_selection,
                                  fitness_source, fitness_migr,
                                  max_fitness_source, max_fitness_migr);
        }
       // Rcout << "mate\n";
        Fish kid = mate(parent1, parent2, size_in_morgan);

       // Rcout << "update fitness and max fitness\n";
        double fit = -2.0;
        if(use_selection) fit = calculate_fitness(kid, select, multiplicative_selection);
        if(fit > new_max_fitness) new_max_fitness = fit;

        new_fitness.push_back(fit);
        new_generation.push_back(kid);
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
                                      bool progress_bar,
                                      arma::mat& frequencies,
                                      bool track_frequency,
                                      const NumericVector& track_markers,
                                      bool track_junctions,
                                      std::vector<double>& junctions,
                                      bool multiplicative_selection,
                                      int num_alleles,
                                      const std::vector<int>& founder_labels,
                                      double migration_rate) {

    // Rcout << "simulate_population: " << multiplicative_selection << "\n";

    bool test = TRUE;

    bool use_selection = FALSE;
    if(select(1, 1) >= 0) use_selection = TRUE;

    double expected_max_fitness = 1e-6;
    std::vector<Fish> pop_1 = source_pop_1;
    std::vector<Fish> pop_2 = source_pop_2;

    double max_fitness_pop_1, max_fitness_pop_2;
    std::vector<double> fitness_pop_1, fitness_pop_2;

    // Rcout << "calculating fitness if necessary\n";

    if(use_selection) {
        for(int j = 0; j < select.nrow(); ++j) {
            if(select(j, 4) < 0) break; // these entries are only for tracking, not for selection calculations
            double local_max_fitness = 0.0;
            for(int i = 1; i < 4; ++i) {
                if(select(j, i) > local_max_fitness) {
                    local_max_fitness = select(j, i);
                }
            }
            expected_max_fitness += local_max_fitness;
        }

        for(auto it = pop_1.begin(); it != pop_1.end(); ++it){
            double fit = calculate_fitness((*it), select, multiplicative_selection);
            if(fit > max_fitness_pop_1) max_fitness_pop_1 = fit;

            if(fit > (expected_max_fitness)) { // little fix to avoid numerical problems
                Rcout << "Expected maximum " << expected_max_fitness << " found " << fit << "\n";
                Rcpp::stop("ERROR in calculating fitness, fitness too large\n");
            }

            fitness_pop_1.push_back(fit);
        }

        for(auto it = pop_2.begin(); it != pop_2.end(); ++it){
            double fit = calculate_fitness((*it), select, multiplicative_selection);
            if(fit > max_fitness_pop_2) max_fitness_pop_2 = fit;

            if(fit > (expected_max_fitness)) { // little fix to avoid numerical problems
                Rcout << "Expected maximum " << expected_max_fitness << " found " << fit << "\n";
                Rcpp::stop("ERROR in calculating fitness, fitness too large\n");
            }

            fitness_pop_2.push_back(fit);
        }
    }

    int updateFreq = total_runtime / 20;
    if(updateFreq < 1) updateFreq = 1;

    if(progress_bar) {
        Rcout << "0--------25--------50--------75--------100\n";
        Rcout << "*";
    }
    R_FlushConsole();

    for(int t = 0; t < total_runtime; ++t) {
      // add code in the future to track junctions
       // if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));
// Rcout << "track frequency is used: " << track_frequency << "\n";
        if(track_frequency) {
            // Rcout << "track frequency start\n"; R_FlushConsole();

            arma::mat local_mat = update_all_frequencies_tibble_dual_pop (pop_1, pop_2,
                                                                          track_markers,
                                                                          founder_labels,
                                                                          t);

            int num_founder_labels = founder_labels.size();
            int num_markers = track_markers.size();
            int local_mat_size = num_founder_labels * num_markers * 2;


            int start_add_marker = local_mat_size * t;

            // Rcout << num_founder_labels << "\t" << num_markers << "\t" << local_mat_size << "\t" << start_add_marker << "\t" << t << "\n";
            // Rcout << "starting appending to large frequencies tibble\n";
            for(int j = 0; j < local_mat_size; ++j) {
              for(int k = 0; k < 5; ++k) {
                // Rcout << t << "\t" << j << "\t" << k << "\t" << start_add_marker + j << "\n";
                // Rcout << frequencies.n_rows << "\t" << frequencies.n_cols << "\n";
                // Rcout << local_mat.n_rows << "\t" << local_mat.n_cols << "\n";
                // Rcout << "source: " << local_mat(j, k) << "\n";
                // Rcout << "target: " << frequencies(start_add_marker + j, k) << "\n";
                frequencies(start_add_marker + j, k)  = local_mat(j, k);
              }
            }
            // Rcout << "done appending\n";
        }

        double new_max_fitness_pop_1, new_max_fitness_pop_2;
        std::vector<double> new_fitness_pop_1, new_fitness_pop_2;
        assert(pop_size.size() == 2);
        // Rcout << "starting nex_pop_migr pop 1\n";
        std::vector<Fish> new_generation_pop_1 = next_pop_migr(pop_1, pop_2, pop_size[0],
                                                               fitness_pop_1, fitness_pop_2,
                                                               max_fitness_pop_1, max_fitness_pop_2,
                                                               select, use_selection, multiplicative_selection,
                                                               migration_rate, new_fitness_pop_1, new_max_fitness_pop_1,
                                                               morgan);
        // Rcout << "starting nex_pop_migr pop 2\n";
        std::vector<Fish> new_generation_pop_2 = next_pop_migr(pop_2, pop_1, pop_size[1],
                                                               fitness_pop_2, fitness_pop_1,
                                                               max_fitness_pop_2, max_fitness_pop_1,
                                                               select, use_selection, multiplicative_selection,
                                                               migration_rate, new_fitness_pop_2, new_max_fitness_pop_2,
                                                               morgan);
        // Rcout << "updating vectors\n";
        pop_1 = new_generation_pop_1;
        pop_2 = new_generation_pop_2;
        fitness_pop_1 = new_fitness_pop_1;
        fitness_pop_2 = new_fitness_pop_2;
        max_fitness_pop_1 = new_max_fitness_pop_1;
        max_fitness_pop_2 = new_max_fitness_pop_2;

        if(t % updateFreq == 0 && progress_bar) {
            Rcout << "**";
        }

        // Rcout << "checking for fixation\n";
        if(is_fixed(pop_1) && is_fixed(pop_2)) {
            Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
            R_FlushConsole();
            std::vector< std::vector < Fish > > output;
            output.push_back(pop_1);
            output.push_back(pop_2);
            return(output);
        }

        Rcpp::checkUserInterrupt();
    }
    if(progress_bar) Rcout << "\n";
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
                            bool progress_bar,
                            bool track_frequency,
                            NumericVector track_markers,
                            bool track_junctions,
                            bool multiplicative_selection,
                            double migration_rate)
{
    // Rcout << "simulate_cpp: " << multiplicative_selection << "\n"; R_FlushConsole();
    bool test = TRUE;

    std::vector< Fish > Pop_1;
    std::vector< Fish > Pop_2;
    int number_of_alleles = -1;
    std::vector<int> founder_labels;

    if(input_population_1[0] > -1e4) {
        Rcout << "Found input populations! converting!\n";  R_FlushConsole();
        Pop_1 = convert_NumericVector_to_fishVector(input_population_1);
        Pop_2 = convert_NumericVector_to_fishVector(input_population_2);

        for(auto it = Pop_1.begin(); it != Pop_1.end(); ++it) {
            update_founder_labels((*it).chromosome1, founder_labels);
            update_founder_labels((*it).chromosome2, founder_labels);
        }

        for(auto it = Pop_2.begin(); it != Pop_2.end(); ++it) {
            update_founder_labels((*it).chromosome1, founder_labels);
            update_founder_labels((*it).chromosome2, founder_labels);
        }

        number_of_alleles = founder_labels.size();
       //  Rcout << "Number of alleles is " << number_of_alleles << "\n";  R_FlushConsole();
    } else {
       //  Rcout << "generating initial population\n";  R_FlushConsole();

        for(int j = 0; j < 2; ++j) {
            for(int i = 0; i < pop_size[j]; ++i) {
                NumericVector focal_freqs = starting_frequencies(j, _);

                int founder_1 = draw_random_founder(focal_freqs);
                int founder_2 = draw_random_founder(focal_freqs);

                Fish p1 = Fish( founder_1 );
                Fish p2 = Fish( founder_2 );

                if(j == 0) Pop_1.push_back(mate(p1,p2, morgan));
                if(j == 1) Pop_2.push_back(mate(p1,p2, morgan));
            }
        }
        for(int i = 0; i < starting_frequencies.ncol(); ++i) {
          founder_labels.push_back(i);
        }
        number_of_alleles = founder_labels.size();
       //  Rcout << "done\n"; R_FlushConsole();
    }

    // Rcout << "Preparing frequencies_table\n"; R_FlushConsole();
    int number_of_markers = track_markers.size();
    // 5 columns: time, loc, anc, type, population
    arma::mat frequencies_table(number_of_markers * number_of_alleles * total_runtime * 2, 5);
    arma::mat initial_frequencies = update_all_frequencies_tibble_dual_pop(Pop_1,
                                                                           Pop_2,
                                                                           track_markers,
                                                                           founder_labels,
                                                                           0);
    // Rcout << "collected initial frequencies\n"; R_FlushConsole();

    std::vector<double> junctions;
    Rcout << "starting simulation\n"; R_FlushConsole();
    std::vector< std::vector< Fish> > output_populations;

    output_populations = simulate_two_populations(Pop_1,
                                                  Pop_2,
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
                                                  migration_rate);
    Rcout << "finished simulation\n";
    arma::mat final_frequencies = update_all_frequencies_tibble_dual_pop(output_populations[0],
                                                               output_populations[1],
                                                               track_markers,
                                                               founder_labels,
                                                               total_runtime);
    // Rcout << "simulation done\n"; R_FlushConsole();
   // if(test) {
   //     return(List::create( Named("population_1")

    //));
   // }

    return List::create( Named("population_1") = convert_to_list(output_populations[0]),
                         Named("population_2") = convert_to_list(output_populations[1]),
                         Named("frequencies") = frequencies_table,
                         Named("initial_frequencies") = initial_frequencies,
                         Named("final_frequencies") = final_frequencies,
                         Named("junctions") = junctions);
}