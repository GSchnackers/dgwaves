
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fluxVector, int uNum){

    std::size_t i, j, k, l;

    std::fill(fluxVector.begin(), fluxVector.end(), 0);

    std::vector<double> scalarProds(uNum * frontierElement.elementTag.size() * frontierElement.numGp, 0);
    std::vector<double> integrations(uNum * frontierElement.nodeTags.size(), 0);

    // Computation of all scalar products at the Gauss points of each frontier element.
    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numGp; ++j)
            for(k = 0; k < uNum; ++k)
                for(l = 0; l < 3; ++l)
                {
                    int scalIndex = i * frontierElement.numGp * uNum + j * uNum + k;
                    int fluxIndex = i * frontierElement.numGp * uNum * 3 + j * uNum * 3 + k * 3 + l;
                    int normIndex = i * frontierElement.numGp * 3 + j * 3 + l;

                    scalarProds[scalIndex] += flux.num[fluxIndex] * frontierElement.normals[normIndex];
                }

    // Computation of the integration vector.
    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
            for(k = 0; k < frontierElement.numGp; ++k)
            {
                int jacobIndex = i * frontierElement.numGp + k;
                int shapeIndex = k * frontierElement.numNodes + j;

                for(l = 0; l < uNum; ++l)
                {
                    int intIndex = i * frontierElement.numNodes * uNum + j * uNum + l;
                    int gpIndex  = i * frontierElement.numGp * uNum + k * uNum + l;
                    

                    integrations[intIndex] += scalarProds[gpIndex] * frontierElement.jacobiansDet[jacobIndex] * \
                                              frontierElement.shapeFunctionsParam[shapeIndex] * \
                                              frontierElement.gaussPointsParam[k * 4 + 3];
                }
            }


    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
        {
            int nodeIndex  = i * frontierElement.numNodes + j;
            
            for(k = 0; k < uNum; ++k)
            {
                int intIndex   = i * frontierElement.numNodes * uNum + j * uNum + k;
                int mainIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes * uNum + \
                                 frontierElement.nodeCorrespondance[nodeIndex].first * uNum + k;

                fluxVector[mainIndex1] += integrations[intIndex];

                if(frontierElement.neighbours[i].second >= 0)
                {
                    int mainIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes * uNum + \
                                     frontierElement.nodeCorrespondance[nodeIndex].second * uNum + k;

                    fluxVector[mainIndex2] -= integrations[intIndex];
                }
            }
        }
            
}