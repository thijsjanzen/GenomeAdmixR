// Copyright 2018 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//

#pragma once

#include <RcppParallel.h>

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
