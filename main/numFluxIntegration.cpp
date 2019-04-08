
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fVector){

    std::size_t i, j, k, l;

    fVector.resize(mainElement.nodeTags.size(), 0);

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numGp; ++j)
            for(k = 0; k < frontierElement.numNodes; ++k)
            {
                int frontNodeIndex = i * frontierElement.numNodes + k;
                int frontGpIndex = i * frontierElement.numGp + j;
                int mainNodeIndex1 = frontierElement.neighbours[frontNodeIndex].first + \
                                    frontierElement.nodeCorrespondance[frontNodeIndex].first; 
                int mainNodeIndex2 = frontierElement.neighbours[frontNodeIndex].second + \
                                    frontierElement.nodeCorrespondance[frontNodeIndex].second;

                for(l = 0; l < 3; ++l)
                {
                    int fluxNormIndex = i * frontierElement.numGp * 3 + j * 3 + l;
                    fVector[mainNodeIndex1] += flux.numGp[fluxNormIndex].first * \
                                               frontierElement.normals[fluxNormIndex];
                    fVector[mainNodeIndex2] += flux.numGp[fluxNormIndex].second * \
                                               frontierElement.normals[fluxNormIndex];
                }
            }
    
}