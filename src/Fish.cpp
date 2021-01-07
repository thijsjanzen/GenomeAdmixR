//
//  Fish.cpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#include "Fish.h"
#include "random_functions.h"
#include <algorithm>


#include <Rcpp.h>
using namespace Rcpp;

std::vector<junction> recombine_new(const std::vector<junction>& chromosome1,
                                    const std::vector<junction>& chromosome2,
                                    const std::vector<double>& recom_positions)
{
    static thread_local auto tl_go = decltype(chromosome1){};
    assert(!chromosome1.genome.empty());    // not strictly enforced by code
    assert(!chromosome2.genome.empty());    // not strictly enforced bu code

    // we need something that is cheaply swappable:
    auto* g1 = &chromosome1;
    auto* g2 = &chromosome2;
    auto& go = tl_go;   // offspring genome: recycle what's already there...
    go.clear();

    // predicate for lower_bound
    auto less = [](const auto& j, double p) { return j.pos < p; };

    // helper lambda to get the value just *before* it.
    // we store the value to the right of a recombination-point but we need the value to the left:
    auto value_at = [](auto begin, auto it) { return (begin != it) ? (it - 1)->right : -1; };

    double left_pos = 0.0;
    auto go_val = -1;
    for (auto right_pos : recom_positions) {
        auto it = std::lower_bound(g1->cbegin(), g1->cend(), left_pos, less);
        auto last = std::lower_bound(it, g1->cend(), right_pos, less);
        // [g1.first, it) : part of the genome *before* left_pos.
        // [it, last) : part of the genome *after or equal to* left_pos but *before* right_pos.
        auto g1_val = value_at(g1->cbegin(), it);
        if (g1_val != go_val) {
            if (it == last || it->pos != left_pos) {
                go.emplace_back(left_pos, g1_val);   // insert change to match
            }
            else {
                ++it;    // corner case: skip spurious double-change
            }
        }
        go.insert(go.end(), it, last);      // append [it, last)
        go_val = value_at(go.begin(), go.end());
        std::swap(g1, g2);
        left_pos = right_pos;
    }
    go.emplace_back(1.0, -1);
    return go;
}

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

void recombine_old(std::vector<junction>& offspring,
                   const std::vector<junction>& chromosome1,
                   const std::vector<junction>& chromosome2,
                   const std::vector<double>& recomPos) {


    std::vector < std::vector<junction>::const_iterator > iters =
        { chromosome1.begin(), chromosome2.begin() };

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
            junction new_junction(recomPos[recompos_cnt], prev_iter->right);
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

    return;
}

std::vector<double> generate_recomPos(int number_of_recombinations,
                                      rnd_t& rndgen) {

    std::vector<double> recomPos(number_of_recombinations, 0);
    for(int i = 0; i < number_of_recombinations; ++i) {
        recomPos[i] = rndgen.uniform();
    }
    std::sort(recomPos.begin(), recomPos.end() );

    if(recomPos.back() < 1.0)
        recomPos.push_back(1.0);

    return recomPos;
}

void Recombine(      std::vector<junction>& offspring,
                     const std::vector<junction>& chromosome1,
                     const std::vector<junction>& chromosome2,
                     double MORGAN,
                     rnd_t& rndgen)  {

    int numRecombinations = rndgen.poisson(MORGAN);

    if (numRecombinations == 0) {
        offspring.insert(offspring.end(),
                         chromosome1.begin(),
                         chromosome1.end());

        return;
    }

    std::vector<double> recomPos = generate_recomPos(numRecombinations,
                                                     rndgen);

#ifdef __unix__
    offspring = recombine_new(chromosome1, chromosome2, recomPos);
#else
    recombine_old(offspring,
                  chromosome1,
                  chromosome2,
                  recomPos);
#endif

    return;
}

Fish mate(const Fish& A, const Fish& B, double numRecombinations,
          rnd_t& rndgen)
{
    Fish offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear(); //just to be sure.

    //first the father chromosome
    int event = rndgen.random_number(2);
    switch(event) {
    case 0:  {
        Recombine(offspring.chromosome1, A.chromosome1, A.chromosome2, numRecombinations, rndgen);
        break;
    }
    case 1: {
        Recombine(offspring.chromosome1, A.chromosome2, A.chromosome1, numRecombinations, rndgen);
        break;
    }
    }


    //then the mother chromosome
    event = rndgen.random_number(2);
    switch(event) {
    case 0:  {
        Recombine(offspring.chromosome2, B.chromosome1, B.chromosome2, numRecombinations, rndgen);
        break;
    }
    case 1: {
        Recombine(offspring.chromosome2, B.chromosome2, B.chromosome1, numRecombinations, rndgen);
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

bool Fish::operator ==(const Fish& other) const {
    if (chromosome1.size() != other.chromosome1.size())
        return false;
    if (chromosome2.size() != other.chromosome2.size())
        return false;

    for (size_t i = 0; i < chromosome1.size(); ++i) {
        if (chromosome1[i] != other.chromosome1[i])
            return false;
    }
    for (size_t i = 0; i < chromosome2.size(); ++i) {
        if (chromosome2[i] != other.chromosome2[i])
            return false;
    }

    return true;
}

bool Fish::operator !=(const Fish& other) const {
    return ! ((*this) == other);
}


Fish::Fish(Fish&& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
}

Fish& Fish::operator=(Fish&& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
    }
    return *this;
}

Fish::Fish(const Fish& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
}

Fish& Fish::operator=(const Fish& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
    }
    return *this;
}