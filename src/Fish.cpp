//
//  Fish.cpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#include "Fish.h"
#include "random_functions.h"
//#include "randomc.h"
#include <algorithm>


#include <Rcpp.h>
using namespace Rcpp;

void add(std::vector< junction>& offspring,
         const junction& new_junction) {

    if (offspring.empty()) {
        offspring.push_back(new_junction);
        return;
    }

    if (new_junction.pos > offspring.back().pos &&
        new_junction.right != offspring.back().right) {
        offspring.push_back(new_junction);
    }
    return;
}

bool do_recombination(std::vector<junction>& offspring,
                      const std::vector<junction>& chromosome1,
                      const std::vector<junction>& chromosome2,
                      std::vector<double>& recomPos) {

    std::vector < std::vector<junction>::const_iterator > iters =
        { chromosome1.begin(), chromosome2.begin() };

    recomPos.push_back(1.0); // for completeness

    int index = 0;
    int recompos_cnt = 0;

    while(true) {

        if ( iters[index]->pos > recomPos[recompos_cnt]  ) {
            // encountered junction point
            // create junction
            index = 1 - index;
            while( iters[index]->pos < recomPos[recompos_cnt]) {
                iters[index]++;
            }

            auto prev_iter = iters[index];
            prev_iter--;
            assert(prev_iter->pos < recomPos[recompos_cnt]);
            junction new_junction(recomPos[recompos_cnt], prev_iter->right);
            // offspring.push_back(new_junction);
            add(offspring, new_junction);

            recompos_cnt++;
        } else {
            add(offspring, (*iters[index]));
            iters[index]++;
        }

        if (offspring.back().right == -1) {
            break;
        }
    }

    return true;
}

std::vector<double> generate_recomPos(int number_of_recombinations) {

    std::vector<double> recomPos(number_of_recombinations, 0);
    for(int i = 0; i < number_of_recombinations; ++i) {
        recomPos[i] = uniform();
    }
    std::sort(recomPos.begin(), recomPos.end() );
    recomPos.erase(std::unique(recomPos.begin(), recomPos.end()), recomPos.end());

    while (recomPos.size() < number_of_recombinations) {
        double pos = uniform();
        recomPos.push_back(pos);
        // sort them, in case they are not sorted yet
        // we need this to remove duplicates, and later
        // to apply crossover
        std::sort(recomPos.begin(), recomPos.end() );
        // remove duplicate recombination sites
        recomPos.erase(std::unique(recomPos.begin(), recomPos.end()), recomPos.end());
    }
    return recomPos;
}

void Recombine(      std::vector<junction>& offspring,
               const std::vector<junction>& chromosome1,
               const std::vector<junction>& chromosome2,
               double MORGAN)  {

    int numRecombinations = poisson_preset();

    if (numRecombinations == 0) {
        offspring.insert(offspring.end(),
                         chromosome1.begin(),
                         chromosome1.end());

        return;
    }

    std::vector<double> recomPos = generate_recomPos(numRecombinations);

    bool recomPos_is_unique = do_recombination(offspring,
                                               chromosome1,
                                               chromosome2,
                                               recomPos);
    // very rarely, the recombination positions are exactly
    // on existing junctions - this should not happen.
    while(recomPos_is_unique == false) {

        recomPos = generate_recomPos(numRecombinations);

        recomPos_is_unique = do_recombination(offspring,
                                              chromosome1,
                                              chromosome2,
                                              recomPos);
    }

    return;
}

Fish mate(const Fish& A, const Fish& B, double numRecombinations)
{
    Fish offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear(); //just to be sure.

    //first the father chromosome
    int event = random_number(2);
    switch(event) {
        case 0:  {
            Recombine(offspring.chromosome1, A.chromosome1, A.chromosome2, numRecombinations);
            break;
        }
        case 1: {
            Recombine(offspring.chromosome1, A.chromosome2, A.chromosome1, numRecombinations);
            break;
        }
    }


    //then the mother chromosome
    event = random_number(2);
    switch(event) {
        case 0:  {
            Recombine(offspring.chromosome2, B.chromosome1, B.chromosome2, numRecombinations);
            break;
        }
        case 1: {
            Recombine(offspring.chromosome2, B.chromosome2, B.chromosome1, numRecombinations);
            break;
        }
    }

    return offspring;
}

Fish::Fish(){

}

junction::junction(){

}

junction::junction(long double loc, int B)  {
    pos = loc;
    right = B;
}

junction::junction(const junction& other) {
    pos = other.pos;
    right = other.right;
}

bool junction::operator ==(const junction& other) const {
    if(pos != other.pos) return false;
    if(right != other.right) return false;

    return true;
}

bool junction::operator <(const junction& other) const {
    return(pos < other.pos);
}

bool junction::operator !=(const junction& other) const {
    return( !( (*this) == other) );
}

Fish::Fish(int initLoc)    {
    junction left(0.0, initLoc);
    junction right(1,  -1);
    chromosome1.push_back( left  );
    chromosome1.push_back( right );
    chromosome2.push_back( left  );
    chromosome2.push_back( right );
}