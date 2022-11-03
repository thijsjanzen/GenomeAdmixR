//
//  Fish.hpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#ifndef Fish_hpp
#define Fish_hpp

#include <stdio.h>
#include <vector>
#include "random_functions.h"

struct junction {
    long double pos;
    int right;

    junction();
    junction(long double loc, int B) ;
    junction(const junction& other);
    bool operator ==(const junction& other) const;
    bool operator !=(const junction& other) const;
};


struct Fish {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish();

    Fish(int initLoc);

    Fish(const Fish& other);
    Fish(Fish&& other);
    Fish& operator=(Fish&& other);
    Fish& operator=(const Fish& other);
};


Fish mate(const Fish& A, const Fish& B, double numRecombinations, rnd_t& rndgen);

#endif /* Fish_hpp */
