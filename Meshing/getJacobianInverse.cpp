#include <cstdio>
#include <iostream>
#include <vector>
#include "meshing.hpp"
#include "structures.hpp"

// Function that inverse all jacobians of elements defined by "element".
void getJacobiansInverse(Element & element){

    std::size_t i, j;

    // Gives the matrix of the invere jacobians the same size as the matrix of jacobians.
    element.jacobiansInverse.resize(element.jacobians.size());

    for(i = 0; i < element.jacobians.size(); i += 9){

        std::vector<double> tmp(9), tmpInverse;

        for(j = 0; j < 9; ++j) tmp[j] = element.jacobians[i + j]; 

        invert(tmp, tmpInverse);

        for(j = 0; j < 9; ++j) element.jacobiansInverse[i + j] = tmpInverse[j];

    }

}