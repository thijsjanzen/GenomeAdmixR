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
#include <algorithm>  // std::sort, std::lower_bound
#include <numeric>    // std::accumulate
#include <stdexcept>  // std::runtime_error
#include <cmath>
#include <cstdint>

// Random number generator. It is always seeded explicitly: the master
// generator of a simulation is seeded from R's RNG (see seed_from_r() in
// util.h), so set.seed() makes simulations reproducible.
struct rnd_t {
  std::mt19937_64 rndgen_;
  std::uniform_real_distribution<> unif_dist =
                                 std::uniform_real_distribution<>(0, 1.0);

  explicit rnd_t(uint64_t seed) : rndgen_(seed) {}

  // independent stream 'stream' derived from 'seed', used to give every
  // block of individuals its own generator in multi-threaded code.
  rnd_t(uint64_t seed, uint64_t stream) {
    std::seed_seq ss{static_cast<uint32_t>(seed),
                     static_cast<uint32_t>(seed >> 32),
                     static_cast<uint32_t>(stream),
                     static_cast<uint32_t>(stream >> 32)};
    rndgen_.seed(ss);
  }

  void set_seed(uint64_t s) {
    rndgen_.seed(s);
  }

  // draw a fresh seed, e.g. to seed the worker generators of one generation
  uint64_t next_seed() {
    return rndgen_();
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
};

// Recombination map for sequence data: entry i is the (relative) probability
// that a crossover falls between marker i - 1 and marker i.
struct emp_genome {
  std::vector< double > cdf_;

 private:
  static void check_weights(const std::vector<double>& weights) {
    if (weights.empty()) {
      throw std::runtime_error("recombination map is empty");
    }
    for (const double weight : weights) {
      if (!(weight >= 0.0) || !std::isfinite(weight)) {
        throw std::runtime_error("recombination map contains negative or "
                                 "non-finite values, are the markers sorted?");
      }
    }
  }

  void build_cdf_from_weights(const std::vector<double>& weights,
                              double total_weight) {
    cdf_.resize(weights.size());
    double cumulative_probability = 0.0;
    for (size_t i = 0; i < weights.size(); ++i) {
      cumulative_probability += weights[i] / total_weight;
      cdf_[i] = cumulative_probability;
    }
  }

  // no information: every interval between markers is equally likely
  void build_uniform_cdf(size_t number_of_markers) {
    cdf_.resize(number_of_markers);
    for (size_t i = 0; i < number_of_markers; ++i) {
      cdf_[i] = number_of_markers > 1 ?
                  static_cast<double>(i) / (number_of_markers - 1) : 1.0;
    }
  }

 public:
  emp_genome() = default;

  explicit emp_genome(const std::vector<double>& weights) {
    check_weights(weights);
    const double total_weight = std::accumulate(weights.begin(),
                                                weights.end(),
                                                0.0);
    if (total_weight > 0.0) {
      build_cdf_from_weights(weights, total_weight);
    } else {
      build_uniform_cdf(weights.size());
    }
    cdf_.back() = 1.0;  // guard against rounding error
  }

  // Default map: crossover probability proportional to the physical
  // distance between adjacent markers. Positions must be sorted.
  static emp_genome from_positions(const std::vector<double>& positions) {
    if (positions.empty()) {
      throw std::runtime_error("positions is empty");
    }
    std::vector<double> distance_to_previous(positions.size(), 0.0);
    for (size_t i = 1; i < positions.size(); ++i) {
      distance_to_previous[i] = positions[i] - positions[i - 1];
    }
    return emp_genome(distance_to_previous);
  }

  // index of the first cdf entry that is not smaller than 'probability'
  size_t index_from_cdf(double probability) const {
    if (cdf_.empty()) {
      throw std::runtime_error("empty cdf");
    }
    auto position = std::lower_bound(cdf_.begin(), cdf_.end(), probability);
    if (position == cdf_.end()) {
      --position;
    }
    return static_cast<size_t>(std::distance(cdf_.begin(), position));
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
