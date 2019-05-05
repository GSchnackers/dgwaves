
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fluxVector){

    std::size_t i, j = 0, k, l;

    std::fill(fluxVector.begin(), fluxVector.end(), 0);

    std::vector<double> scalarProds(6 * frontierElement.elementTag.size() * frontierElement.numGp, 0);
    std::vector<double> integrations(6 * frontierElement.nodeTags.size(), 0);

    // Computation of all scalar products at the Gauss points of each frontier element.
    for(i = 0; i < flux.num.size(); ++i)
        scalarProds[i/3] += flux.num[i] * frontierElement.normals[i/18 + (i % 3)];

    // Computation of the integration vector.
    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
            for(k = 0; k < frontierElement.numGp; ++k)
            {
                int jacobIndex = i * frontierElement.numGp + k;
                int shapeIndex = k * frontierElement.numNodes + j;

                for(l = 0; l < 6; ++l)
                {
                    int intIndex = i * frontierElement.numNodes * 6 + j * 6 + l;
                    int gpIndex = i * frontierElement.numGp * 6 + k * 6 + l;
                    

                    integrations[intIndex] += scalarProds[gpIndex] * frontierElement.jacobiansDet[jacobIndex] * \
                                              frontierElement.shapeFunctionsParam[shapeIndex] * \
                                              frontierElement.gaussPointsParam[k * 4 + 3];
                }
            }


    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
            for(k = 0; k < 6; ++k)
            {
                int intIndex   = i * frontierElement.numNodes * 6 + j * 6 + k;
                int nodeIndex  = i * frontierElement.numNodes + j;
                int mainIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes * 6 + \
                                 frontierElement.nodeCorrespondance[nodeIndex].first * 6 + k;

                fluxVector[mainIndex1] += integrations[intIndex];

                if(frontierElement.neighbours[i].second >= 0)
                {
                    int mainIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes * 6 + \
                                     frontierElement.nodeCorrespondance[nodeIndex].second * 6 + k;

                    fluxVector[mainIndex2] -= integrations[intIndex];
                }
            }
            
}