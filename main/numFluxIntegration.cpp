
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fluxVector){

    std::size_t i, j, k, l;

    std::fill(fluxVector.begin(), fluxVector.end(), 0);

    std::vector<double> scalarProds(frontierElement.elementTag.size() * frontierElement.numGp, 0);
    std::vector<double> integrations(frontierElement.elementTag.size() * frontierElement.numNodes, 0);

    // Computation of all scalar products at the Gauss points of each frontier element.
    for(i = 0; i < flux.num.size(); ++i)
    {
        scalarProds[i/3] += flux.num[i] * frontierElement.normals[i];
    }

    // Computation of the integration vector.
    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
        {
            int index = i * frontierElement.numNodes + j;

            for(k = 0; k < frontierElement.numGp; ++k)
            {
                int shapeIndex = k * frontierElement.numNodes + j;
                int gpIndex = i * frontierElement.numGp + k;

                integrations[index] += scalarProds[gpIndex]  * frontierElement.jacobiansDet[gpIndex] *\
                                       frontierElement.shapeFunctionsParam[shapeIndex] * \
                                       frontierElement.gaussPointsParam[k * 4 + 3];
            }

        }

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
        {
            int index = i * frontierElement.numNodes + j;
            int mainIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes + \
                             frontierElement.nodeCorrespondance[index].first;

            fluxVector[mainIndex1] += integrations[index];

            if(frontierElement.neighbours[i].second >= 0)
            {
                int mainIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes + \
                                 frontierElement.nodeCorrespondance[index].second;

                fluxVector[mainIndex2] -= integrations[index];
            }
        }
            
}