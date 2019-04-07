
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fMatrix){

    std::size_t i, j, k;

    fMatrix.resize(mainElement.nodeTags.size());

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numGp; ++j)
            for(k = 0; k < mainElement.numNodes; ++k)
            {
                int indexF1 = frontierElement.neighbours[i].first * mainElement.numNodes;
                int indexF2 = frontierElement.neighbours[i].second * mainElement.numNodes;
                int indexNode1 = frontierElement.neighbours[i].first * mainElement.numNodes + k;
                int indexNode2 = frontierElement.neighbours[i].second * mainElement.numNodes + k;



            }
}