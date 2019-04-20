#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void timeMarching(const Element & mainElement, const std::vector<double> & SFProd, \
                  const std::vector<double> & fluxVector, const double step, const double t,
                  std::vector<double> & kVector){

    std::size_t i, j, k;
    std::fill(kVector.begin(), kVector.end(), 0);

    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
        {
            int uIndex = i * mainElement.numNodes + j;
            double tmpProd = 0;

            for(k = 0; k < mainElement.numNodes; ++k)
            {
        
                int matrixIndex = i * mainElement.numNodes * mainElement.numNodes + \
                                  j * mainElement.numNodes + k;

                int vecIndex = i * mainElement.numNodes + k;

                kVector[uIndex] += mainElement.massMatrixInverse[matrixIndex] * (SFProd[vecIndex] - fluxVector[vecIndex]);
                
            }

        } 

}