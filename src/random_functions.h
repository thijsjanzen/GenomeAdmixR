//
//  random_functions.hpp
//
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#ifndef random_functions_hpp
#define random_functions_hpp

#include <random>
#include <vector>

double uniform();
int random_number(int n);
void set_seed(unsigned seed);

int poisson_preset();
void set_poisson(double lambda);

void fill_cum_marker_dist(const std::vector<double>& positions);
std::vector< size_t > recompos();

int num_mutations();
void set_mutation_rate(double mu,
                       int num_markers);


#endif /* random_functions_hpp */
