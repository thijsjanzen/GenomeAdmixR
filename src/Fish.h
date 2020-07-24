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

struct junction {
    long double pos;
    int right;

    junction();
    junction(long double loc, int B) ;
    junction(const junction& other);
    bool operator ==(const junction& other) const;
    bool operator <(const junction& other) const;
    bool operator !=(const junction& other) const;
};


struct Fish {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish();

    Fish(int initLoc);
};


Fish mate(const Fish& A, const Fish& B, double numRecombinations);

#endif /* Fish_hpp */
