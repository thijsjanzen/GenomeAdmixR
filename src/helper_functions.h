//
//  helper_functions.hpp
//
//
//  Created by Thijs Janzen on 20/11/2018.
//

#ifndef helper_functions_hpp
#define helper_functions_hpp

#include <vector>
#include "Fish.h"
#include "random_functions.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;


// bool verify_individual_cpp(const Fish& Nemo);
// bool verify_pop_cpp(const std::vector< Fish >& pop);
bool matching_chromosomes(const std::vector< junction >& v1,
                          const std::vector< junction >& v2);

template <typename FISH>
bool is_fixed(const std::vector< FISH >& v);

template <typename FISH>
NumericVector update_frequency(const std::vector< FISH >& v,
                               double m,
                               int num_alleles);

template <typename FISH>
arma::mat update_all_frequencies(const std::vector< FISH >& pop,
                                 const NumericVector& markers,
                                 int number_of_founders);

template <typename FISH>
double calc_mean_junctions(const std::vector< FISH> & pop);

int draw_prop_fitness(const std::vector<double>& fitness,
                      double maxFitness);

template <typename FISH>
std::vector< FISH > convert_NumericVector_to_fishVector(const NumericVector& v);

List convert_to_list(const std::vector<Fish<chromosome_junctions>>& v);

template <typename FISH>
double calculate_fitness(const FISH& focal,
                         const NumericMatrix& select,
                         bool multiplicative_selection);

int draw_random_founder(const NumericVector& v);
void update_founder_labels(const std::vector<junction> chrom,
                           std::vector<int>& founder_labels);


template <typename FISH>
arma::mat update_frequency_tibble(const std::vector< FISH >& v,
                                  double m,
                                  const std::vector<int>& founder_labels,
                                  int t,
                                  double morgan);

template <typename FISH>
arma::mat update_all_frequencies_tibble(const std::vector< FISH >& pop,
                                        const NumericVector& markers,
                                        const std::vector<int>& founder_labels,
                                        int t,
                                        double morgan);

template <typename FISH>
arma::mat update_all_frequencies_tibble_dual_pop(const std::vector< FISH >& pop_1,
                                                 const std::vector< FISH >& pop_2,
                                                 const NumericVector& markers,
                                                 const std::vector<int>& founder_labels,
                                                 int t,
                                                 double morgan);

template <typename FISH>
arma::mat update_frequency_tibble_dual_pop(const std::vector< FISH >& pop_1,
                                           const std::vector< FISH >& pop_2,
                                           double marker,
                                           const std::vector<int>& founder_labels,
                                           int t,
                                           double morgan);

NumericVector scale_markers(const Rcpp::NumericVector& markers,
                            double morgan);

#endif /* helper_functions_hpp */
