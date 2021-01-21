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

double calc_mean_junctions(const std::vector< Fish<chromosome_junctions> > & pop);

int draw_prop_fitness(const std::vector<double>& fitness,
                      double maxFitness);

std::vector< Fish<chromosome_junctions> >
  convert_NumericVector_to_fishVector(const Rcpp::NumericVector& v);

List convert_to_list(const std::vector<Fish<chromosome_junctions>>& v);

int draw_random_founder(const NumericVector& v);
void update_founder_labels(const chromosome_junctions& chrom,
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

double calc_fitness(const Fish<chromosome_junctions> & focal,
                    const Rcpp::NumericMatrix& select,
                    bool multiplicative_selection);

#endif /* helper_functions_hpp */
