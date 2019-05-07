#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void timeMarching(const Element & mainElement, const std::vector<double> & SFProd, \
                  const std::vector<double> & fluxVector, std::vector<double> & kVector, int uNum){

    std::size_t i, j, k, l;
    std::fill(kVector.begin(), kVector.end(), 0);

    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < mainElement.numNodes; ++k)
                for(l = 0; l < uNum; ++l)
                {
                    int uIndex      = i * mainElement.numNodes * uNum + j * uNum + l;
                    int matrixIndex = i * mainElement.numNodes * mainElement.numNodes + \
                                      j * mainElement.numNodes + k;

                    int vecIndex = i * mainElement.numNodes * uNum + k * uNum + l;

                    kVector[uIndex] += mainElement.massMatrixInverse[matrixIndex] * \
                                       (SFProd[vecIndex] - fluxVector[vecIndex]);
                    
                }
            


}