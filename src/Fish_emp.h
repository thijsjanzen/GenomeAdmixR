#ifndef Fish_emp_hpp
#define Fish_emp_hpp

#include <stdio.h>
#include <vector>
#include "random_functions.h"
#include <RcppArmadillo.h>

struct Fish_emp {
  std::vector< int > chromosome1;
  std::vector< int > chromosome2;

  Fish_emp() {
  }

  Fish_emp(const std::vector< int >& c1,
           const std::vector< int >& c2) :
    chromosome1(c1), chromosome2(c2) {
  }
};

inline std::vector< int > gamete(double morgan,
                          rnd_t& rndgen,
                          const emp_genome& emp_gen,
                          const Fish_emp& parent)  {

  std::vector<int> recom_pos = emp_gen.recompos(morgan,
                                                rndgen);

  std::vector< int > recombined_chromosome;

  if (recom_pos.size() <= 1) {
    recombined_chromosome = rndgen.random_number(2) ? parent.chromosome1 : parent.chromosome2;
    return recombined_chromosome;
  }

  int index = rndgen.random_number(2);
  int prev_start = 0;

  for(size_t i = 0; i < recom_pos.size(); ++i) {
    auto start = prev_start;
    auto end   = recom_pos[i];

    for (int j = start; j < end; ++j) {
      // this code can be more optimized, but first we go for something that works.
      if (j < 0 || j > parent.chromosome1.size()) {
        Rcpp::stop("recombine out of range");
      }
      if (index == 0) recombined_chromosome.push_back(parent.chromosome1[j]);
      if (index == 1) recombined_chromosome.push_back(parent.chromosome2[j]);
    }

    prev_start = end;
    index = 1 - index;
  }

  return recombined_chromosome;
}

#endif