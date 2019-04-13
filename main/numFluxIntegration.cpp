
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fluxVector){

    std::size_t i, j, k, l;

    std::vector<double> intTemp(frontierElement.nodeTags.size(), 0); // variable that stores for each node the evaluation of the integral.

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
        {

            int nodeIndex = i * frontierElement.numNodes + j;

            int vecIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes + \
                                    frontierElement.nodeCorrespondance[nodeIndex].first;
            int vecIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes + \
                                    frontierElement.nodeCorrespondance[nodeIndex].second;


            for(k = 0; k < frontierElement.numGp; ++k)
            {
                
                int shapeIndex = k * frontierElement.numNodes + j;

                int detIndex = i * frontierElement.numGp + k;

                double prodScal = 0;

                for(l = 0; l < 3; ++l)
                {
                    
                    int fluxIndex = i * frontierElement.numGp * 3 + k * 3 + l;

                    prodScal += flux.num[fluxIndex] * frontierElement.normals[fluxIndex];

                }
                    
                intTemp[nodeIndex] += prodScal * \
                                       frontierElement.shapeFunctionsParam[shapeIndex] * \
                                       frontierElement.gaussPointsParam[4 * k + 3] * 
                                       frontierElement.jacobiansDet[detIndex];
 
            }

            fluxVector[vecIndex1] += intTemp[nodeIndex];
            if(frontierElement.neighbours[i].second >= 0) fluxVector[vecIndex2] -= intTemp[nodeIndex];

        }

        
}