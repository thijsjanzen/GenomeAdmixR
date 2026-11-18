//  Copyright (c) 2026, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

// bare minimal tbb-stub header.
// just enough to to let *this* package pass the
// RCPP_PARALLEL_USE_TBB shenanigan on Alpine Linux

#pragma once

#include <cstdlib>
#include <Rcpp.h>
#include <RcppParallel.h>   // pull RCPP_PARALLEL_USE_TBB
#include <utility>


#if RCPP_PARALLEL_USE_TBB == 0

// everything looks so single-threaded here :(

#include <algorithm>

namespace tbb {

class task_arena {
public:
  explicit task_arena(int /*num_threads*/) {}
  template<typename Func> void execute(const Func& f) { f(); }
  static const int automatic = -1;
};


template<typename T>
class blocked_range {
  T begin_;
  T end_;
public:
  blocked_range(T begin, T end) : begin_(begin), end_(end) {}
  T begin() const { return begin_; }
  T end() const { return end_; }
};


template<typename T, typename Func>
inline void parallel_for(blocked_range<T> b, const Func& f) {
  f(b);
}

}  // namespace tbb

#endif  // RCPP_PARALLEL_USE_TBB == 0


#include <algorithm>
#include <cstdint>
#include "random_functions.h"

// Draw a 64 bit seed from R's random number generator, such that set.seed()
// in R makes simulations reproducible. The Rcpp generated wrappers in
// RcppExports.cpp hold an RNGScope, so R's RNG state is loaded here.
// Written by Claude (Anthropic).
inline uint64_t seed_from_r() {
  const double two_to_the_32 = 4294967296.0;
  const auto high_bits = static_cast<uint64_t>(R::unif_rand() * two_to_the_32);
  const auto low_bits  = static_cast<uint64_t>(R::unif_rand() * two_to_the_32);
  return (high_bits << 32) ^ low_bits;
}

// Calls task(index, random_generator) for every index in [0, number_of_tasks),
// possibly in parallel.
//
// The indices are split into fixed blocks of 64, and every block gets its own
// random number generator, derived from (seed, block number). The outcome
// therefore depends only on 'seed', not on the number of threads or on the
// order in which threads run.
//
// num_threads < 1 means: use all available cores.
// 'task' must not call the R API; report errors with C++ exceptions instead.
// Written by Claude (Anthropic).
template <typename Task>
void parallel_blocks(size_t number_of_tasks,
                     int num_threads,
                     uint64_t seed,
                     const Task& task) {
  constexpr size_t block_size = 64;
  const size_t number_of_blocks = (number_of_tasks + block_size - 1) / block_size;
  const int maximum_concurrency = num_threads < 1 ? -1 : num_threads;  // -1: automatic

  auto run_block = [&](size_t block_number) {
    rnd_t random_generator(seed, block_number);
    const size_t first_index = block_number * block_size;
    const size_t end_index = std::min(number_of_tasks, first_index + block_size);
    for (size_t index = first_index; index < end_index; ++index) {
      task(index, random_generator);
    }
  };

  auto run_range_of_blocks = [&](const tbb::blocked_range<size_t>& block_range) {
    for (size_t block_number = block_range.begin();
         block_number < block_range.end();
         ++block_number) {
      run_block(block_number);
    }
  };

  tbb::task_arena arena(maximum_concurrency);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, number_of_blocks),
                      run_range_of_blocks);
  });
}
