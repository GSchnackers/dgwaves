
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fVector){

    std::size_t i, j, k, l;

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
            for(k = 0; k < frontierElement.numGp; ++k)
            {
                int gaussIndex = i * frontierElement.numGp + k;
                int shapeIndex = k * frontierElement.numGp + j;
                
                int dVecIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes + \
                                    frontierElement.nodeCorrespondance[gaussIndex].first;
                int dVecIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes + \
                                    frontierElement.nodeCorrespondance[gaussIndex].second;

                fVector[dVecIndex1] = fVector[dVecIndex2] =  0;

                for(l = 0; l < 3; ++l)
                {
                    
                    int fluxIndex = i * frontierElement.numGp * 3 + k * 3 + l;
                    
                    fVector[dVecIndex1] += flux.numGp[fluxIndex].first * \
                                           frontierElement.normals[fluxIndex] * \
                                           frontierElement.shapeFunctionsParam[shapeIndex] * \
                                           frontierElement.gaussPointsParam[4 * k + 3];

                    if(frontierElement.neighbours[i].second >= 0)
                        fVector[dVecIndex2] += flux.numGp[fluxIndex].second * \
                                               frontierElement.normals[fluxIndex] * \
                                               frontierElement.shapeFunctionsParam[shapeIndex] * \
                                               frontierElement.gaussPointsParam[4 * k + 3];

                }
            }

        
}