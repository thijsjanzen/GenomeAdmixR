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
#endif /* random_functions_hpp */
