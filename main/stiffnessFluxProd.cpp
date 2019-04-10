#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// This function computes the product of the stiffness matrices and the physical flux vector at nodal values.
void stiffnessFluxProd(const Element & mainElement, const Quantity & flux, std::vector<double> & prod){

    std::size_t i, j, k;
    
    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
        {
            int prodIndex = i * mainElement.numNodes + j;
            prod[prodIndex] = 0;

            for(k = 0; k < mainElement.numNodes; ++k)
            {
                int stiffIndex = i * mainElement.numNodes * mainElement.numNodes + k * mainElement.numNodes + \
                                 + k;

                int vecIndex = i * mainElement.numNodes * 3 + k;
                

                prod[prodIndex] += mainElement.stiffnessMatrixX[stiffIndex] * flux.node[vecIndex] + \
                                   mainElement.stiffnessMatrixY[stiffIndex] * flux.node[vecIndex + 1] + \
                                   mainElement.stiffnessMatrixZ[stiffIndex] * flux.node[vecIndex + 2];
            }
        }

}