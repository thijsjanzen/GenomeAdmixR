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

bool do_recombination(std::vector<junction>& offspring,
                      const std::vector<junction>& chromosome1,
                      const std::vector<junction>& chromosome2,
                      const std::vector<double>& recomPos) {

    std::vector< junction > toAdd; //first create junctions on exactly the recombination positions
    for(int i = 0; i < recomPos.size(); ++i) {
        junction temp;
        temp.right = -1.0;
        temp.pos = recomPos[i];
        toAdd.push_back(temp);
    }

    for(auto i = (chromosome1.begin()+1); i != chromosome1.end(); ++i) {
        long double leftpos = (*(i-1)).pos;
        long double rightpos = (*i).pos;

        for(int j = 0; j < recomPos.size(); ++j) {
            if(recomPos[j] == leftpos) {
                return false;
            }
            if(recomPos[j] == rightpos) {
                return false;
            }

            if(recomPos[j] > leftpos) {
                if(recomPos[j] < rightpos) {
                    if(j % 2 == 0) { // even, so chrom1 = L, chrom2 = R
                        if(j < toAdd.size()) {
                            //       toAdd[j].left = (*i).left;
                        }
                    } else { // uneven so chrom1 = R, chrom2 = L
                        if(j < toAdd.size()) {
                            toAdd[j].right = (*(i-1)).right;
                        }
                    }
                }
            }
        }
    }

    for(auto i = (chromosome2.begin()+1); i != chromosome2.end(); ++i) {
        long double leftpos = (*(i-1)).pos;
        long double rightpos = (*i).pos;

        for(int j = 0; j < recomPos.size(); ++j) {
            if(recomPos[j] == leftpos) {
                return false;
            }
            if(recomPos[j] == rightpos) {
                return false;
            }

            if(recomPos[j] > leftpos) {
                if(recomPos[j] < rightpos) {
                    if(j % 2 == 0) { //even, so chrom1 = L, chrom2 = R
                        if(j < toAdd.size()) {
                            toAdd[j].right = (*(i-1)).right;
                        }
                    } else { //uneven so chrom1 = R, chrom2 = L
                        if(j < toAdd.size()) {
                            //    toAdd[j].left = (*i).left;
                        }
                    }
                }
            }
        }
    }

    for(int i = 0; i < toAdd.size(); ++i) {
        if(toAdd[i].right == -1 && toAdd[i].pos < 1) {
            Rcout << "This break point was not addressed!\n";
            stop("Error in toAdd");
        }
        offspring.push_back(toAdd[i]);
    }

    //now we have to add the other junctions from chrom1 and chrom2.
    long double leftpos = 0;
    long double rightpos = 0;


    for(int i = 0; i < (recomPos.size() + 1); ++i) {
        rightpos = 1.0;
        if(i < recomPos.size()) rightpos = recomPos[i];

        if(i % 2 == 0) { //even, so take from chromosome 1
            for(auto it = chromosome1.begin(); it != chromosome1.end(); ++it) {
                if((*it).pos >= leftpos && (*it).pos <= rightpos) {
                    offspring.push_back((*it));
                }
                if((*it).pos > rightpos) { //we are past the recombination section
                    break;
                }
            }
        }

        if(i % 2 == 1) { //odd, so take from chromosome 2
            for(auto it = chromosome2.begin(); it != chromosome2.end(); ++it) {
                if((*it).pos >= leftpos && (*it).pos <= rightpos) {
                    offspring.push_back((*it));
                }
                if((*it).pos > rightpos) { //we are past the recombination section
                    break;
                }
            }
        }

        //move forward
        leftpos = rightpos;
    }

    std::sort(offspring.begin(), offspring.end());
    // gatekeeper code to not allow false junctions to be introduced

    std::vector<junction> temp_offspring = offspring;
    offspring.clear();
    for(int i = 0; i < temp_offspring.size(); ++i) {  // extra checks to make sure no memory access errors
        bool add = true;

        if(i > 0) {
            if(temp_offspring[i].right == temp_offspring[i-1].right) add = false;

            if(temp_offspring[i].pos == temp_offspring[i-1].pos) add = false;
        }

        // if a bad memory access call happens, we get weird numbers in .right:
        if(abs(temp_offspring[i].right) > 1000) add = false;

        if(add == true) {
            offspring.push_back(temp_offspring[i]);
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