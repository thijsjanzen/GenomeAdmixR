//
//  random_functions.hpp
//
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#pragma once

#include <random>
#include <vector>
#include <algorithm> // std::sort
#include <thread>
#include <chrono>
#include <iostream>
#include <cmath>


struct rnd_t {
  std::mt19937_64 rndgen_;
  std::uniform_real_distribution<> unif_dist = std::uniform_real_distribution<>(0, 1.0);
  std::uniform_int_distribution<> rand_num_dist;

  rnd_t() {
    auto seed = get_seed();
    rndgen_ = std::mt19937_64(seed);
  }

  rnd_t(unsigned int seed) {
    auto local_seed = get_seed() + seed;
    rndgen_ = std::mt19937_64(local_seed);
  }

  void set_seed(unsigned int s) {
    rndgen_ = std::mt19937_64(s);
  }

  double uniform() {
    return unif_dist(rndgen_);
  }

  int random_number(int n) {
    if (n < 1) return 0;
    return std::uniform_int_distribution<int> (0, n - 1)(rndgen_);
  }

  size_t poisson(double lambda) {
    return std::poisson_distribution<size_t>(lambda)(rndgen_);
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

  emp_genome() {
  }

  emp_genome(const emp_genome& other) {
    cdf_ = other.cdf_;
  }

  emp_genome& operator=(const emp_genome& other) {
    if (this != &other) {
      cdf_ = other.cdf_;
    }
    return *this;
  }

  template <typename T>
  emp_genome(const std::vector<T>& positions) {
    if (positions.empty()) {
      throw std::runtime_error("positions is empty");
    }
    auto total_sum = std::accumulate(positions.begin(),
                                positions.end(),
                                0.0);

    double s = 0.0;
    double mult = 1.0 / total_sum;
    cdf_.resize(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
      s += positions[i] * mult;
      cdf_[i] = s;
    }
    return;
  }

  size_t index_from_cdf(double p) const {

    if (cdf_.empty()) throw "empty cdf";
    
    if (cdf_.back() <= 0.0) return static_cast<size_t>(p * cdf_.size());
    // find index belonging to p
    return static_cast<size_t>(std::distance(cdf_.begin(),
                                             std::lower_bound(cdf_.begin(),
                                                              cdf_.end(),
                                                              p)));
  }

  std::vector< size_t > recompos(double morgan,
                                 rnd_t& rndgen) const {
    size_t num_break_points = rndgen.poisson(morgan);
    std::vector< size_t > indices;
    for(size_t i = 0; i < num_break_points; ++i) {
      auto found_index = index_from_cdf(rndgen.uniform());
      if (found_index > 0) { // first position can not be recombination point
        indices.push_back(found_index);
      }
    }
    std::sort(indices.begin(), indices.end());
    indices.push_back(cdf_.size());
    return indices;
  }
};
