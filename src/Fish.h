//
//  Fish.hpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#ifndef Fish_hpp
#define Fish_hpp

#include <vector>
#include "random_functions.h"

// A junction marks the start of an ancestry block: from 'pos' onwards
// (up to the next junction, or the end of the chromosome at 1.0) the
// ancestry is 'right'. Every chromosome starts with a junction at 0.0;
// there is no end sentinel.
struct junction {
    double pos = 0.0;
    int right = -1;

    junction() = default;
    junction(double loc, int B) : pos(loc), right(B) {}

    bool operator ==(const junction& other) const {
      return pos == other.pos && right == other.right;
    }
    bool operator !=(const junction& other) const {
      return !(*this == other);
    }
};


struct Fish {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish() = default;
    explicit Fish(int initLoc);
};


Fish mate(const Fish& A, const Fish& B, double morgan, rnd_t& rndgen);

#endif /* Fish_hpp */
