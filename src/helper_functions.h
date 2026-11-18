//
//  helper_functions.hpp
//
//
//  Created by Thijs Janzen on 20/11/2018.
//

#ifndef helper_functions_hpp
#define helper_functions_hpp

#include <vector>
#include <array>
#include "Fish.h"
#include "Fish_emp.h"
#include "random_functions.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void force_output();

int draw_prop_fitness(const std::vector<double>& fitness,
                      const double& maxFitness,
                      rnd_t& rndgen);

std::vector< std::array<double, 5> > convert_select_from_r(const Rcpp::NumericMatrix& select,
                                                           double scale = 1.0);

//// ancestry (Fish) functions ////

int get_ancestry(const std::vector< junction >& chrom, double marker);

int draw_random_founder(const NumericVector& v,
                        rnd_t& rndgen);

void update_founder_labels(const std::vector<junction> chrom,
                           std::vector<int>& founder_labels);

std::vector< Fish > convert_NumericVector_to_fishVector(const NumericVector& v);

List convert_to_list(const std::vector<Fish>& v);

double calculate_fitness(const Fish& focal,
                         const std::vector<std::array<double, 5>>& select,
                         bool multiplicative_selection);

double total_junctions(const std::vector< Fish >& pop);

arma::mat update_frequency_tibble(const std::vector< Fish >& v,
                                  double m,
                                  const std::vector<int>& founder_labels,
                                  int t,
                                  double morgan);

arma::mat update_all_frequencies_tibble(const std::vector< Fish >& pop,
                                        const std::vector<double>& markers,
                                        const std::vector<int>& founder_labels,
                                        int t,
                                        double morgan);

NumericVector scale_markers(const Rcpp::NumericVector& markers,
                            double morgan);


//// sequence (Fish_emp) functions ////

std::vector< Fish_emp > convert_numeric_matrix_to_fish_vector(
    const Rcpp::NumericMatrix& input_population) ;

double calculate_fitness(const Fish_emp& focal,
                         const std::vector<std::array<double, 5>>& select,
                         const std::vector<double>& locations,
                         bool multiplicative_selection);

std::vector< std::vector<double > > update_frequency_tibble(const std::vector< Fish_emp >& pop,
                                                            size_t marker_index,
                                                            double pos,
                                                            int t);

arma::mat update_all_frequencies_tibble(const std::vector< Fish_emp >& pop,
                                        const std::vector<double>& markers,
                                        const std::vector<double>& locations,
                                        int t,
                                        double morgan);

List convert_to_list(const std::vector<Fish_emp>& v,
                     const std::vector<double>& locations);

int find_location(const std::vector<double>& markers,
                  double pos);

void mutate(Fish_emp& indiv,
            const std::vector<std::vector<double>>& sub_matrix,
            const double& mutation_rate,
            rnd_t& rndgen);

#endif /* helper_functions_hpp */
