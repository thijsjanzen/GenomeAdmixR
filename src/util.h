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

namespace task_arena {

constexpr size_t automatic = size_t(-1);

}  // namespace task_arena


class global_control {
public:
  enum parameter {
    max_allowed_parallelism,
    thread_stack_size,
    terminate_on_exception
  };

  global_control(parameter /*p*/, size_t /*value*/) {}
  ~global_control() {}
  static size_t active_value(parameter /*param*/);  // undefined
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
inline void parallel_for(blocked_range<T> b, const Func f) {
  f(b);
}

}  // namespace tbb


// function name is lying.
inline size_t get_rcpp_num_threads() {
  return 1;
}

// do nothing in this case.
inline void set_num_threads() {
  return;
}




#else  // if RCPP_PARLLEL_USE_TBB = 0


// probably the cleanest way to retrieve RcppParallel's concurrency setting
// set by RcppParallel::setThreadOptions(numThreads)
inline size_t get_rcpp_num_threads() {
  auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  return (nullptr == nt_env)
    ? tbb::task_arena::automatic  // -1
  : static_cast<size_t>(std::atoi(nt_env));
}

inline void set_num_threads() {
  auto num_threads = get_rcpp_num_threads();
  auto global_control =
    tbb::global_control(tbb::global_control::max_allowed_parallelism,
                        num_threads);
}

#endif