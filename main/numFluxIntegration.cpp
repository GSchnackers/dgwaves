
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
    {
        if(!(i % 18) && i) j += 3;
        scalarProds[i/3] += flux.num[i] * frontierElement.normals[j + (i % 3)];
    }

    // Computation of the integration vector.
    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < 6; ++j)
            for(k = 0; k < frontierElement.numNodes; ++k)
            {
                int index = i * frontierElement.numNodes * 6 + j + k * 6;

                for(l = 0; l < frontierElement.numGp; ++l)
                {

                    int shapeIndex = l * frontierElement.numNodes + k;
                    int gpIndex = i * frontierElement.numGp * 6 + j + l * 6;
                    int jacobIndex = i * frontierElement.numGp + l;

                    integrations[index] += scalarProds[gpIndex] * frontierElement.jacobiansDet[jacobIndex] * \
                                        frontierElement.shapeFunctionsParam[shapeIndex] * \
                                        frontierElement.gaussPointsParam[l * 4 + 3];
                }

            }

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < 6; ++j)
            for(k = 0; k < frontierElement.numNodes; ++k)
            {
                int index = i * frontierElement.numNodes + k;
                int mainIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes * 6 + j + \
                                frontierElement.nodeCorrespondance[index].first * 6;

                fluxVector[mainIndex1] += integrations[index];

                if(frontierElement.neighbours[i].second >= 0)
                {
                    int mainIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes * 6 + j + \
                                    frontierElement.nodeCorrespondance[index].second * 6;

                    fluxVector[mainIndex2] -= integrations[index];
                }
            }
            
}