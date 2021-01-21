//
//  Fish.cpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#include "Fish.h"
#include "random_functions.h"
#include <algorithm>


#include <Rcpp.h>
using namespace Rcpp;



std::vector<junction> recombine_new(const chromosome_junctions& chromosome1,
                                    const chromosome_junctions& chromosome2,
                                    const std::vector<double>& recom_positions)
{

  static auto tl_go = decltype(chromosome1.genome){};

  // we need something that is cheaply swappable:
  auto* g1 = &chromosome1.genome;
  auto* g2 = &chromosome2.genome;
  auto& go = tl_go;   // offspring genome: recycle what's already there...
  go.clear();

  // predicate for lower_bound
  auto less = [](const auto& j, double p) { return j.pos < p; };

  // helper lambda to get the value just *before* it.
  // we store the value to the right of a recombination-point but we need the value to the left:
  auto value_at = [](auto begin, auto it) { return (begin != it) ? (it - 1)->right : -1; };

  double left_pos = 0.0;
  auto go_val = -1;
  for (auto right_pos : recom_positions) {
    auto it = std::lower_bound(g1->cbegin(), g1->cend(), left_pos, less);
    auto last = std::lower_bound(it, g1->cend(), right_pos, less);
    // [g1.first, it) : part of the genome *before* left_pos.
    // [it, last) : part of the genome *after or equal to* left_pos but *before* right_pos.
    auto g1_val = value_at(g1->cbegin(), it);
    if (g1_val != go_val) {
      if (it == last || it->pos != left_pos) {
        go.emplace_back(left_pos, g1_val);   // insert change to match
      }
      else {
        ++it;    // corner case: skip spurious double-change
      }
    }
    go.insert(go.end(), it, last);      // append [it, last)
    go_val = value_at(go.begin(), go.end());
    std::swap(g1, g2);
    left_pos = right_pos;
  }
  go.emplace_back(1.0, -1);
  return go;
}


std::vector<double> generate_recomPos(size_t number_of_recombinations) {

  std::vector<double> recomPos(number_of_recombinations, 0);
  for (size_t i = 0; i < number_of_recombinations; ++i) {
    recomPos[i] = uniform();
  }
  std::sort(recomPos.begin(), recomPos.end() );
  recomPos.erase(std::unique(recomPos.begin(), recomPos.end()), recomPos.end());

  recomPos.push_back(1.0);
  return recomPos;
}

chromosome_junctions::chromosome_junctions(const chromosome_junctions& chromosome1,
                                           const chromosome_junctions& chromosome2,
                                           double morgan) {

  int num_positions = poisson_preset();
  std::vector < double > recom_positions = generate_recomPos(num_positions);

  // do recombination
  if (uniform() < 0.5) {
    genome = recombine_new(chromosome1,
                         chromosome2,
                         recom_positions);
  } else {
    genome = recombine_new(chromosome2,
                           chromosome1,
                           recom_positions);
  }
  return;
}
