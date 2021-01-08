//
//  random_functions.cpp
//
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#include "random_functions.h"
#include <Rcpp.h>
#include <random>
using namespace Rcpp;

std::random_device rd;
std::mt19937 rndgen(rd());  //< The one and only random number generator
std::uniform_real_distribution<> unif_dist = std::uniform_real_distribution<>(0, 1.0);
std::poisson_distribution<int> poisson_preset_dist;
std::uniform_int_distribution<> rand_num_dist;

int random_number(int n)    {
    return std::uniform_int_distribution<> (0, n-1)(rndgen);
}

double uniform()    {
    return  unif_dist(rndgen);  //std::uniform_real_distribution<>(0, 1.0)(rndgen);
}

int poisson_preset() {
    return poisson_preset_dist(rndgen);
}

void set_poisson(double lambda) {
    poisson_preset_dist = std::poisson_distribution<int>(lambda);
}

void set_seed(unsigned seed)    {
    std::mt19937 new_randomizer(seed);
    rndgen = new_randomizer;
}