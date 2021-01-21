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

struct junction {
    long double pos;
    int right;

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

    chromosome_junctions(const chromosome_junctions& A,
                         const chromosome_junctions& B,
                         double morgan);

    chromosome_junctions& operator=(const chromosome_junctions& other) {
        genome = other.genome;
        return *this;
    }

    chromosome_junctions(const chromosome_junctions& other) {
        genome = other.genome;
    }
};

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
