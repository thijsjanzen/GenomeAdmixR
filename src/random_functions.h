//
//  random_functions.hpp
//
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#ifndef random_functions_hpp
#define random_functions_hpp

#include <random>
#include <vector>

struct rnd_t {
  std::mt19937 rndgen_;

  rnd_t() {
    std::random_device rd;
    rndgen_ = std::mt19937(rd());
  }

  rnd_t(unsigned int seed) {
    rndgen_ = std::mt19937(seed);
  }

  void set_seed(unsigned int s) {
    rndgen_ = std::mt19937(s);
  }

  std::uniform_real_distribution<> unif_dist =
    std::uniform_real_distribution<>(0, 1.0);

  double uniform() {
    return unif_dist(rndgen_);
  }

  int random_number(unsigned int n) {
    return std::uniform_int_distribution<> (0, n-1)(rndgen_);
  }

  int poisson(double lambda) {
    return std::poisson_distribution<int>(lambda)(rndgen_);
  }

  std::vector<double> generate_random_markers(int number_of_markers) {
    if(number_of_markers < 0) {
      // regularly spaced markers
      std::vector<double> markers;
      double di = 1.0 / (1 + (-1 * number_of_markers));
      double marker_pos = di;
      while(marker_pos < 1.f) {
        markers.push_back(marker_pos);
        marker_pos += di;
      }
      return(markers);
    }

    std::vector<double> markers(number_of_markers);
    for(int i = 0; i < number_of_markers; ++i) {
      markers[i] = uniform();
    }
    std::sort(markers.begin(), markers.end());

    markers.erase( unique( markers.begin(), markers.end() ), markers.end() );
    while(markers.size() < number_of_markers) {
      int diff = number_of_markers - markers.size();
      for(int i = 0; i < diff; ++i) {
        markers.push_back(uniform());
      }
      std::sort(markers.begin(), markers.end());
      markers.erase( unique( markers.begin(), markers.end() ), markers.end() );
    }
    return(markers);
  }
};

#endif /* random_functions_hpp */
