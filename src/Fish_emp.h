#ifndef Fish_emp_hpp
#define Fish_emp_hpp

#include <stdio.h>
#include <vector>
#include "random_functions.h"

struct Fish_emp {
  std::vector< int > chromosome1;
  std::vector< int > chromosome2;

  Fish_emp()
    {}

  Fish_emp(const std::vector< int >& c1,
           const std::vector< int >& c2) :
    chromosome1(c1), chromosome2(c2) {
  }

  Fish_emp(Fish_emp&& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
  }

  Fish_emp& operator=(Fish_emp&& other) {
    if (this != &other) {
      chromosome1 = other.chromosome1;
      chromosome2 = other.chromosome2;
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

  std::vector< int > gamete(const std::vector< size_t >& recom_pos) {

    if (recom_pos.size() == 1) {
      if(random_number(2)) {
        return chromosome1;
      }
      return chromosome2;
    }

    std::vector < std::vector<int>::const_iterator > iters = {chromosome1.begin(),
                                                              chromosome2.begin()};
    std::vector< int > recombined_chromosome;
    int index = random_number(2);
    size_t prev_start = 0;

    for(int i = 0; i < recom_pos.size(); ++i) {
      auto start = iters[index] + prev_start;
      auto end   = iters[index] + recom_pos[i];

      prev_start = recom_pos[i];
      recombined_chromosome.insert(recombined_chromosome.end(), start, end);
      index = 1 - index;
    }
    assert(!recombined_chromosome.empty());
    assert(recombined_chromosome.size() > 0);

    return recombined_chromosome;
  }
};

#endif