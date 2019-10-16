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

int poisson(double lambda)    {
    return std::poisson_distribution<int>(lambda)(rndgen);
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

std::vector<double> generate_random_markers(int number_of_markers) {
    if(number_of_markers < 0) {
        // regularly spaced markers
        std::vector<double> markers;
        double di = 1.0 / (1 + (-1 * number_of_markers));
        double marker_pos = di;
        while(marker_pos < 1.f) {
            markers.push_back(marker_pos);
            marker_pos += di;
        }
        return(markers);
    }

    std::vector<double> markers(number_of_markers);
    for(int i = 0; i < number_of_markers; ++i) {
        markers[i] = uniform();
    }
    std::sort(markers.begin(), markers.end());

    markers.erase( unique( markers.begin(), markers.end() ), markers.end() );
    while(markers.size() < number_of_markers) {
        int diff = number_of_markers - markers.size();
        for(int i = 0; i < diff; ++i) {
            markers.push_back(uniform());
        }
        std::sort(markers.begin(), markers.end());
        markers.erase( unique( markers.begin(), markers.end() ), markers.end() );
    }
    return(markers);
}


