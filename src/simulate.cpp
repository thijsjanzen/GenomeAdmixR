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
                                        const std::vector<int>& founder_labels) {

    //Rcout << "simulate_population: " << multiplicative_selection << "\n";

    bool use_selection = FALSE;
    if(select(1, 1) >= 0) use_selection = TRUE;

    double expected_max_fitness = 1e-6;
    std::vector<Fish> Pop = sourcePop;
    std::vector<double> fitness;
    double maxFitness = -1;

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

        for(auto it = Pop.begin(); it != Pop.end(); ++it){
            double fit = calculate_fitness((*it), select, multiplicative_selection);
            if(fit > maxFitness) maxFitness = fit;

            if(fit > (expected_max_fitness)) { // little fix to avoid numerical problems
                Rcout << "Expected maximum " << expected_max_fitness << " found " << fit << "\n";
                Rcpp::stop("ERROR in calculating fitness, fitness too large\n");
            }

            fitness.push_back(fit);
        }
    }

    int updateFreq = total_runtime / 20;
    if(updateFreq < 1) updateFreq = 1;

    if(progress_bar) {
        Rcout << "0--------25--------50--------75--------100\n";
        Rcout << "*";
    }

    for(int t = 0; t < total_runtime; ++t) {

        if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));

        if(track_frequency) {
            //Rcout << "track frequency start\n";
            for(int i = 0; i < track_markers.size(); ++i) {
                if(track_markers[i] < 0) break;
                arma::mat local_mat = update_frequency_tibble(Pop,
                                                              track_markers[i],
                                                              founder_labels,
                                                              t);

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

        std::vector<Fish> newGeneration;
        std::vector<double> newFitness;
        double newMaxFitness = -1.0;
        //Rcout << "updating fish\n";
        for(int i = 0; i < pop_size; ++i)  {
            int index1 = 0;
            int index2 = 0;
            if(use_selection) {
                index1 = draw_prop_fitness(fitness, maxFitness);
                index2 = draw_prop_fitness(fitness, maxFitness);
                while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness);
            } else {
                index1 = random_number( (int)Pop.size() );
                index2 = random_number( (int)Pop.size() );
                while(index2 == index1) index2 = random_number( (int)Pop.size() );
            }

            Fish kid = mate(Pop[index1], Pop[index2], morgan);

            newGeneration.push_back(kid);

            double fit = -2.0;
            if(use_selection) fit = calculate_fitness(kid, select, multiplicative_selection);
            if(fit > newMaxFitness) newMaxFitness = fit;

            if(fit > expected_max_fitness) {
                Rcout << "Expected maximum " << expected_max_fitness << " found " << fit << "\n";
                Rcpp::stop("ERROR in calculating fitness, fitness too large\n");
            }

            newFitness.push_back(fit);
        }

        if(t % updateFreq == 0 && progress_bar) {
            Rcout << "**";
        }

        if(is_fixed(Pop)) {
            Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
            R_FlushConsole();
            return(Pop);
        }

        Rcpp::checkUserInterrupt();

        Pop = newGeneration;
        newGeneration.clear();
        fitness = newFitness;
        maxFitness = newMaxFitness;
   //     Rcout << "done updating, again!\n";
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
              bool multiplicative_selection)
{
    //Rcout << "simulate_cpp: " << multiplicative_selection << "\n";

    std::vector< Fish > Pop;
    int number_of_alleles = number_of_founders;
    std::vector<int> founder_labels;

    if(input_population[0] > -1e4) {
        //Rcout << "Found input population! converting!\n";
        Pop = convert_NumericVector_to_fishVector(input_population);

        number_of_founders = 0;
        for(auto it = Pop.begin(); it != Pop.end(); ++it) {
            update_founder_labels((*it).chromosome1, founder_labels);
            update_founder_labels((*it).chromosome2, founder_labels);
        }
        number_of_alleles = founder_labels.size();
        //Rcout << "Number of alleles is " << number_of_alleles << "\n";
    } else {
        std::vector<double> starting_freqs = as< std::vector<double> >(starting_proportions);
        for(int i = 0; i < pop_size; ++i) {
            int founder_1 = draw_random_founder(starting_freqs);
            int founder_2 = draw_random_founder(starting_freqs);

            Fish p1 = Fish( founder_1 );
            Fish p2 = Fish( founder_2 );

            Pop.push_back(mate(p1,p2, morgan));
        }
        for(int i = 0; i < number_of_alleles; ++i) {
            founder_labels.push_back(i);
        }
    }

    arma::mat frequencies_table;

    if(track_frequency) {
        //Rcout << "Preparing frequencies_table\n";
       // int number_entries = track_markers.size();
        int number_of_markers = track_markers.size();
        //arma::cube x(total_runtime, number_of_alleles, number_entries); // n_row, n_col, n_slices, type
        arma::mat x(number_of_markers * number_of_alleles * total_runtime, 4); // 4 columns: time, loc, anc, type
        frequencies_table = x;
    }

    arma::mat initial_frequencies = update_all_frequencies_tibble(Pop, track_markers, founder_labels, 0);
  //  Rcout << "collected initial frequencies\n";

    std::vector<double> junctions;
  //  Rcout << "starting simulation\n";
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
                                                      founder_labels);
    //Rcout << "finished simulation\n";
    arma::mat final_frequencies = update_all_frequencies_tibble(outputPop, track_markers, founder_labels, total_runtime);

    return List::create( Named("population") = convert_to_list(outputPop),
                         Named("frequencies") = frequencies_table,
                         Named("initial_frequencies") = initial_frequencies,
                         Named("final_frequencies") = final_frequencies,
                         Named("junctions") = junctions);
}

// [[Rcpp::export]]
void test_fish_functions() {
    Fish test_fish;

    junction temp;
    junction temp2(0.5, 0);
    junction temp3(0.5, 0);
    if(temp2 == temp3) {
        temp = temp2;
    }

    bool temp400 = (temp2 != temp3);
    Rcout << temp400 << "\t" << "this is only for testing\n";

    junction temp4(temp);

    test_fish.chromosome1.push_back(temp);
    test_fish.chromosome1.push_back(temp2);
    test_fish.chromosome1.push_back(temp3);
    test_fish.chromosome1.push_back(temp4);

    Fish test_fish2 = test_fish;

    if(test_fish == test_fish2) {
        Rcout << "fishes are equal!\n";
    }

    Fish test_fish3(5);
    bool b = (test_fish == test_fish3);
    Rcout << b << "\t" << "this is only for testing\n";

    std::vector< junction > chrom;
    chrom.push_back(temp);

    Fish test_fish4(chrom, chrom);

    std::vector< Fish > pop;
    pop.push_back(test_fish);
    pop.push_back(test_fish2);
    pop.push_back(test_fish3);
    pop.push_back(test_fish4);

    verify_pop_cpp(pop);


    return;
}
