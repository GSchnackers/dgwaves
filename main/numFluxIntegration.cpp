
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fluxVector){

    std::size_t i, j, k, l;

    std::vector<double> intVector(frontierElement.nodeTags.size(), 0);

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
        {
            int frontNodeIndex = i * frontierElement.numNodes + j;
            int mainIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes +\
                             frontierElement.nodeCorrespondance[frontNodeIndex].first;
            int mainIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes +\
                             frontierElement.nodeCorrespondance[frontNodeIndex].second;
            std::cout << mainIndex1 << " " << mainIndex2 << std::endl;

            for(k = 0; k < frontierElement.numGp; ++k)
            {
                double prodScal = 0;
                int shapeIndex = k * frontierElement.numNodes + j;
                int detIndex = i * frontierElement.numGp + k;

                for(l = 0; l < 3; ++l)
                {
                    int vecGpIndex = i * frontierElement.numGp * 3 + k * 3 + l;
                    prodScal += flux.num[vecGpIndex] * frontierElement.normals[vecGpIndex];
                }

                 intVector[frontNodeIndex] += prodScal * frontierElement.shapeFunctionsParam[shapeIndex] * \
                                              frontierElement.gaussPointsParam[k * 4 + 3] * \
                                              frontierElement.jacobiansDet[detIndex];


            }

            fluxVector[mainIndex1] += intVector[frontNodeIndex];

            if(frontierElement.neighbours[i].second >= 0)
            {
               fluxVector[mainIndex2] -= intVector[frontNodeIndex];
            }
            


        }
            

        
}