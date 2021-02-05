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
#include <algorithm> // std::sort
#include <thread>
#include <chrono>
#include <iostream>
#include <cmath>

struct rnd_t {
  std::mt19937 rndgen_;
  std::uniform_real_distribution<> unif_dist = std::uniform_real_distribution<>(0, 1.0);
  std::poisson_distribution<int> poisson_preset_dist;
  std::vector< double > cum_marker_dist;
  std::uniform_int_distribution<> rand_num_dist;
  std::binomial_distribution<int> mutate_num;

  rnd_t() {
    auto seed = get_seed();
  //  std::cerr << "initializing rnd_t with: " << seed << "\n" << std::flush;
    rndgen_ = std::mt19937(seed);
  }

  rnd_t(unsigned int seed) {
    rndgen_ = std::mt19937(seed);
  }

  void set_seed(unsigned int s) {
    rndgen_ = std::mt19937(s);
  }

  double uniform() {
    return unif_dist(rndgen_);
  }

  int random_number(unsigned int n) {
    return std::uniform_int_distribution<> (0, n-1)(rndgen_);
  }

  int poisson(double lambda) {
    return std::poisson_distribution<int>(lambda)(rndgen_);
  }

  void set_poisson(double lambda) {
    poisson_preset_dist = std::poisson_distribution<int>(lambda);
  }

  int poisson_preset() {
    return poisson_preset_dist(rndgen_);
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
    return indices;
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

  void set_mutation_rate(double mu,
                         int num_markers) {
    mutate_num = std::binomial_distribution<int>(num_markers, mu);
  }

  int num_mutations() {
    return mutate_num(rndgen_);
  }

  int get_seed() {
    const auto tt = static_cast<int64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    auto tid = std::this_thread::get_id();
    const uint64_t e3{ std::hash<std::remove_const_t<decltype(tid)>>()(tid) };
    auto output = static_cast<int>(tt + e3);
    if (output < 0) output *= -1;
    return output;
  }
};

#endif /* random_functions_hpp */
