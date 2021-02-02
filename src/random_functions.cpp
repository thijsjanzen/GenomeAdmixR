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
std::vector< double > cum_marker_dist;

std::binomial_distribution<int> mutate_num;


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


void fill_cum_marker_dist(const std::vector<double>& positions) {
    auto total_sum = std::accumulate(positions.begin(),
                                     positions.end(), 0.0);
    double s = 0.0;
    double mult = 1.0 / total_sum;
    cum_marker_dist.clear();
    cum_marker_dist.resize(positions.size());
    for (int i = 0; i < positions.size(); ++i) {
        s += positions[i] * mult;
        cum_marker_dist[i] = s;
    }
    return;
}

size_t index_from_cdf(double p) {
    // find index belonging to p
    return static_cast<size_t>(std::distance(cum_marker_dist.begin(),
                                             std::lower_bound(cum_marker_dist.begin(),
                                                              cum_marker_dist.end(),
                                                              p)));
}

std::vector< size_t > recompos() {
    int num_break_points = poisson_preset();
    std::vector< size_t > indices;
    for(size_t i = 0; i < num_break_points; ++i) {
        auto found_index = index_from_cdf(uniform());
        if (found_index > 0)
            indices.push_back(found_index);
    }
    std::sort(indices.begin(), indices.end());
    indices.push_back(cum_marker_dist.size());
  //  for(size_t i = 0; i < indices.size(); ++i) {
//        Rcout << indices[i] << " ";
  //  }
//    Rcout << "\n";
    return indices;
}

void set_mutation_rate(double mu,
                       int num_markers) {
  mutate_num = std::binomial_distribution<int>(num_markers, mu);
}

int num_mutations() {
  return mutate_num(rndgen);
}
