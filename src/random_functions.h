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
  std::uniform_int_distribution<> rand_num_dist;

  rnd_t() {
    auto seed = get_seed();
  //  std::cerr << "initializing rnd_t with: " << seed << "\n" << std::flush;
    rndgen_ = std::mt19937(seed);
  }

  rnd_t(unsigned int seed) {
    auto local_seed = get_seed() + seed;
    rndgen_ = std::mt19937(local_seed);
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

  int get_seed() {
    const auto tt = static_cast<int64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    auto tid = std::this_thread::get_id();
    const uint64_t e3{ std::hash<std::remove_const_t<decltype(tid)>>()(tid) };
    auto output = static_cast<int>(tt + e3);
    if (output < 0) output *= -1;
    return output;
  }
};

struct emp_genome {

  std::vector< double > cdf_;
  int total_sum;

  emp_genome() {

  }

  template <typename T>
  emp_genome(const std::vector<T>& positions) {
    total_sum = std::accumulate(positions.begin(),
                                positions.end(), 0.0);
  //  std::cerr << total_sum << "\n";
    double s = 0.0;
    double mult = 1.0 / total_sum;
    cdf_.resize(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
      cdf_[i] = s;
      s += positions[i] * mult;
 //     std::cerr << cdf_[i] << " ";
    }
 //   std::cerr << "\n";
    return;
  }

  std::vector< size_t > recompos(double morgan,
                                 rnd_t& rndgen) const {
    int num_break_points = rndgen.poisson(morgan);
    std::vector< size_t > indices;
    for(size_t i = 0; i < num_break_points; ++i) {
      auto found_index = index_from_cdf(rndgen.uniform());
      if (found_index > 0)
        indices.push_back(found_index);
    }
    std::sort(indices.begin(), indices.end());
    indices.push_back(cdf_.size());
    return indices;
  }


  size_t index_from_cdf(double p) const {
    // find index belonging to p
    return static_cast<size_t>(std::distance(cdf_.begin(),
                                             std::lower_bound(cdf_.begin(),
                                                              cdf_.end(),
                                                              p)));
  }
};

#endif /* random_functions_hpp */
