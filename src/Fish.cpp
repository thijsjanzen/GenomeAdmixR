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
#include <vector>


std::vector<junction> recombine_new(const std::vector<junction>& chromosome1,
                                    const std::vector<junction>& chromosome2,
                                    const std::vector<double>& recom_positions) {
  // we need something that is cheaply swappable:
  auto* g1 = &chromosome1;
  auto* g2 = &chromosome2;
  std::vector<junction> go;

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
  return go;
}

std::vector<double> generate_recomPos(size_t number_of_recombinations,
                                      rnd_t& rndgen) {

  std::vector<double> recomPos(number_of_recombinations, 0);
  for (size_t i = 0; i < number_of_recombinations; ++i) {
    recomPos[i] = rndgen.uniform();
  }
  std::sort(recomPos.begin(), recomPos.end() );
  recomPos.erase(std::unique(recomPos.begin(), recomPos.end()), recomPos.end());

  recomPos.push_back(1.0);
  return recomPos;
}

void Recombine(      std::vector<junction>& offspring,
                     const std::vector<junction>& chromosome1,
                     const std::vector<junction>& chromosome2,
                     double MORGAN,
                     rnd_t& rndgen)  {

  size_t numRecombinations = rndgen.poisson(MORGAN);

  if (numRecombinations == 0) {
    offspring = chromosome1;
    return;
  }

  std::vector<double> recomPos = generate_recomPos(numRecombinations, rndgen);

  offspring = recombine_new(chromosome1, chromosome2, recomPos);

  return;
}

Fish mate(const Fish& A, const Fish& B, double morgan,
          rnd_t& rndgen) {
  Fish offspring;

  // gamete from parent A
  if (rndgen.random_number(2) == 0) {
    Recombine(offspring.chromosome1, A.chromosome1, A.chromosome2, morgan, rndgen);
  } else {
    Recombine(offspring.chromosome1, A.chromosome2, A.chromosome1, morgan, rndgen);
  }

  // gamete from parent B
  if (rndgen.random_number(2) == 0) {
    Recombine(offspring.chromosome2, B.chromosome1, B.chromosome2, morgan, rndgen);
  } else {
    Recombine(offspring.chromosome2, B.chromosome2, B.chromosome1, morgan, rndgen);
  }

  return offspring;
}

Fish::Fish(int initLoc) {
  junction left(0.0, initLoc);
  chromosome1.push_back(left);
  chromosome2.push_back(left);
}
