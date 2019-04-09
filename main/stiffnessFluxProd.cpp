#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// This function computes the product of the stiffness matrix and the physical flux vector at nodal values.
void stiffnessFluxProd(const Element & mainElement, const Quantity & flux, std::vector<double> & prod){

    std::size_t i, j, k;

    prod.resize(mainElement.nodeTags.size(), 0);
    
    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < mainElement.numNodes; ++k)
            {
                int stiffIndex = i * mainElement.numNodes * mainElement.numNodes + k * mainElement.numNodes + \
                                 + k;

                int vecIndex = i * mainElement.numNodes * 3 + k;
                int prodIndex = i * mainElement.numNodes + j;

                prod[prodIndex] += mainElement.stiffnessMatrixX[stiffIndex] * flux.node[vecIndex] + \
                                   mainElement.stiffnessMatrixY[stiffIndex] * flux.node[vecIndex + 1] + \
                                   mainElement.stiffnessMatrixZ[stiffIndex] * flux.node[vecIndex + 2];
            }

}