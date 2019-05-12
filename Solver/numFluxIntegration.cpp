
#include <cstdio>
#include <iostream>
#include "structures.hpp"
#include <OMP.h>

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fluxVector, int uNum){

    std::size_t i, j, k, l;

    std::fill(fluxVector.begin(), fluxVector.end(), 0);

    // Computation of the integration vector.
    #pragma omp parallel for shared(i, flux, frontierElement, fluxVector)
    for(i = 0; i < frontierElement.elementTag.size(); ++i)
    {
        int iGp = i * frontierElement.numGp;
        int iGpU = iGp * uNum;
        int mainElemIdx1 = frontierElement.neighbours[i].first * mainElement.numNodes * uNum;
        int mainElemIdx2 = frontierElement.neighbours[i].second * mainElement.numNodes * uNum;

        for(j = 0; j < frontierElement.numNodes; ++j)
        { 
            int mainNodeIdx1 = mainElemIdx1 + \
                               frontierElement.nodeCorrespondance[i * frontierElement.numNodes + j].first * uNum;
            int mainNodeIdx2 = mainElemIdx2 + \
                               frontierElement.nodeCorrespondance[i * frontierElement.numNodes + j].second * uNum;

            for(k = 0; k < frontierElement.numGp; ++k)
            {
                int jacobIndex = iGp + k;
                int shapeIndex = k * frontierElement.numNodes + j;

                for(l = 0; l < uNum; ++l)
                {
                    int gpIndex  = iGpU + k * uNum + l;

                    double tmp = flux.num[gpIndex] * frontierElement.jacobiansDet[jacobIndex] * \
                                 frontierElement.shapeFunctionsParam[shapeIndex] * \
                                 frontierElement.gaussPointsParam[k * 4 + 3];

                    fluxVector[mainNodeIdx1 + l] += tmp;

                    if(frontierElement.neighbours[i].second >= 0) fluxVector[mainNodeIdx2 + l] -= tmp;
                }
            }
        }
    }

}