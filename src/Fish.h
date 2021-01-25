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
#include <random>

template <typename CHROMO>
struct Fish {
    CHROMO chromosome1;
    CHROMO chromosome2;

    Fish();

    Fish(const CHROMO& c1, const CHROMO& c2) :
        chromosome1(c1), chromosome2(c2) {
    }

    Fish(int initLoc) :
        chromosome1(CHROMO(initLoc)),
        chromosome2(CHROMO(initLoc)) {
    }

    Fish(const Fish& other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
    }

    Fish(Fish&& other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
    }
    Fish& operator=(Fish&& other) {
        if (this != &other) {
            chromosome1 = other.chromosome1;
            chromosome2 = other.chromosome2;
        }
        return *this;
    }
    Fish& operator=(const Fish& other) {
        if (this != &other) {
            chromosome1 = other.chromosome1;
            chromosome2 = other.chromosome2;
        }
        return *this;
    }

    CHROMO gamete(double morgan) {
        return CHROMO(chromosome1, chromosome2, morgan);
    }
};

#endif /* Fish_hpp */
