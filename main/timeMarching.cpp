#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void timeMarching(const Element & mainElement, const std::vector<double> & SFProd, \
                  const std::vector<double> & fluxVector, std::vector<double> & kVector){

    std::size_t i, j, k, l;
    std::fill(kVector.begin(), kVector.end(), 0);

    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < 6; ++k)
                for(l = 0; l < mainElement.numNodes; ++l)
                {
                    int uIndex = i * mainElement.numNodes * 6 + k + j * 6;
                    int matrixIndex = i * mainElement.numNodes * mainElement.numNodes + \
                                    j * mainElement.numNodes + l;

                    int vecIndex = i * mainElement.numNodes * 6 + k + l * 6;

                    kVector[uIndex] += mainElement.massMatrixInverse[matrixIndex] * \
                                       (SFProd[vecIndex] - fluxVector[vecIndex]);
                    
                }
            


}