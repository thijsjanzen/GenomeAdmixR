#ifndef Fish_emp_hpp
#define Fish_emp_hpp

#include <stdio.h>
#include <vector>
#include "random_functions.h"

struct Fish_emp {
  std::vector< int > chromosome1;
  std::vector< int > chromosome2;

  Fish_emp() {
      chromosome1 = std::vector<int>(); // empty stuff
      chromosome2 = std::vector<int>();
  }

  Fish_emp(Fish_emp&& other) {
    chromosome1 = std::move(other.chromosome1);
    chromosome2 = std::move(other.chromosome2);
  }

  Fish_emp& operator=(Fish_emp&& other) {
    if (this != &other) {
      chromosome1 = std::move(other.chromosome1);
      chromosome2 = std::move(other.chromosome2);
    }
    return *this;
  }

  Fish_emp(const Fish_emp& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
  }

  Fish_emp& operator=(const Fish_emp& other) {
    if (this != &other) {
      chromosome1 = other.chromosome1;
      chromosome2 = other.chromosome2;
    }
    return *this;
  }

  Fish_emp(const std::vector< int >& c1,
           const std::vector< int >& c2) :
    chromosome1(c1), chromosome2(c2) {
  }

  std::vector< int > gamete(double morgan,
                            rnd_t& rndgen,
                            const emp_genome& emp_gen) const {

    std::vector<size_t> recom_pos = emp_gen.recompos(morgan,
                                                     rndgen);

    if (recom_pos.size() <= 1) {
      if(rndgen.random_number(2)) {
        return chromosome1;
      }
      return chromosome2;
    }

    std::vector< int > recombined_chromosome;
    int index = rndgen.random_number(2);
    size_t prev_start = 0;
  //  assert(chromosome1.size() == chromosome2.size());

    for(size_t i = 0; i < recom_pos.size(); ++i) {
      auto start = prev_start;
      auto end   = recom_pos[i];
      prev_start = recom_pos[i];

      if (index == 0) {
        for (size_t j = start; j < end; ++j) {
          recombined_chromosome.push_back(chromosome1[j]);
        }
      } else {
        for (size_t j = start; j < end; ++j) {
          recombined_chromosome.push_back(chromosome2[j]);
        }
      }

      index = 1 - index;
    }

    return recombined_chromosome;
  }
};

#endif