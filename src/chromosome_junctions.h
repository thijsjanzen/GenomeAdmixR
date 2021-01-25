#ifndef chromosome_junctions_hpp
#define chromosome_junctions_hpp

#include "random_functions.h"
#include <vector>

struct junction {
  junction()
  {}

  junction(long double loc, int B)  {
    pos = loc;
    right = B;
  }

  junction(const junction& other) {
    pos = other.pos;
    right = other.right;
  }

  int get_ancestry() const {
    return right;
  }

  bool operator ==(const junction& other) const {
    if(pos != other.pos) return false;
    if(right != other.right) return false;

    return true;
  }

  bool operator <(const junction& other) const {
    return(pos < other.pos);
  }

  bool operator !=(const junction& other) const {
    return( !( (*this) == other) );
  }

// private:  // TODO: make private
  long double pos;
  int right;
};

struct chromosome_junctions {
  std::vector< junction > genome;

  chromosome_junctions(int init) {
    genome.push_back(junction(0.0, init));
    genome.push_back(junction(1.0, -1));
  }

  chromosome_junctions() {
    genome.push_back(junction(0.0, 0.0));
    genome.push_back(junction(1.0, -1));
  }

  chromosome_junctions& operator=(const chromosome_junctions& other) {
    genome = other.genome;
    return *this;
  }

  chromosome_junctions(const chromosome_junctions& other) {
    genome = other.genome;
  }

  int get_num_junctions() const {
    return genome.size() - 2;
  }

  int get_ancestry(double pos) const {
    auto less = [](const auto& j, double p) { return j.pos < p; };

    auto it = std::lower_bound(genome.cbegin(), genome.cend(), pos, less);
    return (it - 1)->right;
  }

  std::vector<junction> recombine(const chromosome_junctions& chromosome1,
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

  chromosome_junctions(const chromosome_junctions& chromosome1,
                       const chromosome_junctions& chromosome2,
                       double morgan) {

    int num_positions = poisson_preset();
    std::vector < double > recom_positions = generate_recomPos(num_positions);

    // do recombination
    if (uniform() < 0.5) {
      genome = recombine(chromosome1,
                         chromosome2,
                         recom_positions);
    } else {
      genome = recombine(chromosome2,
                         chromosome1,
                         recom_positions);
    }
    return;
  }
};


#endif