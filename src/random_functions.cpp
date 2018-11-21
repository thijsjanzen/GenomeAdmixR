//
//  random_functions.cpp
//  
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#include "random_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

double uniform()
{
    return R::runif(0.0, 1.0);
}

long double long_uniform() {
    long double base = 1e9;
    long double a = R::runif(0, 1) * base;
    long double b = R::runif(0, 1);
    long double output = (a+b) / base;
    return(output);
}

int random_number(int n)
{
    return (int)(R::runif(0.0, 1.0 * n));
}

double poisson(double lambda)
{
    return R::rpois(lambda);
}

