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
inline uint64_t seed_from_r() {
  const auto hi = static_cast<uint64_t>(R::unif_rand() * 4294967296.0);
  const auto lo = static_cast<uint64_t>(R::unif_rand() * 4294967296.0);
  return (hi << 32) ^ lo;
}

// Apply body(i, rng) to every i in [0, n), possibly in parallel.
// Work is split into fixed blocks, and every block gets its own random number
// generator derived from (seed, block index). The outcome therefore depends
// only on 'seed', and not on the number of threads or on scheduling.
// num_threads < 1 means: use all available cores.
// 'body' must not call the R API; report errors with C++ exceptions.
template <typename Body>
void parallel_blocks(size_t n, int num_threads, uint64_t seed,
                     const Body& body) {
  constexpr size_t block_size = 64;
  const size_t num_blocks = (n + block_size - 1) / block_size;
  const int concurrency = num_threads < 1 ? -1 : num_threads;

  tbb::task_arena arena(concurrency);
  arena.execute([&] {
    tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_blocks),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t b = r.begin(); b < r.end(); ++b) {
          rnd_t rndgen(seed, b);
          const size_t end = std::min(n, (b + 1) * block_size);
          for (size_t i = b * block_size; i < end; ++i) {
            body(i, rndgen);
          }
        }
      });
  });
}
