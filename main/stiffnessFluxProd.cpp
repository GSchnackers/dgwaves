#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// This function computes the product of the stiffness matrix and the physical flux vector at nodal values.
void stiffnessFluxProd(const Element & mainElement, const Quantity & flux, std::vector<double> & prod){

    std::size_t i, j, k;

    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < mainElement.numNodes; ++k)
            {
                int stiffIndex = i * mainElement.numNodes * mainElement.numNodes + \
                                    j * mainElement.numNodes + k;
                int fluxIndex = i * mainElement.numNodes * 3 + k * 3;
                int prodIndex = i * mainElement.numNodes + j;

                prod[prodIndex] += mainElement.stiffnessMatrixX[stiffIndex] * flux.node[fluxIndex] + \
                                    mainElement.stiffnessMatrixY[stiffIndex] * flux.node[fluxIndex + 1] + \
                                    mainElement.stiffnessMatrixZ[stiffIndex] * flux.node[fluxIndex + 2];

            }

}