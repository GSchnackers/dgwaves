#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "structures.hpp"
#include <omp.h>

// This function computes the product of the stiffness matrices and the physical flux vector at nodal values.
void stiffnessFluxProd(const Element & mainElement, const Quantity & flux, std::vector<double> & prod, int uNum){

    std::size_t i, j, k, l;

    // Components along x, y and z of the gradient.
    std::fill(prod.begin(), prod.end(), 0);
    
    #pragma omp parallel for shared(i, prod, mainElement) private(j,k)
    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < mainElement.numNodes; ++k)
            {
                int stiffIndex = i * mainElement.numNodes * mainElement.numNodes + j * mainElement.numNodes + k;

                for(l = 0; l < uNum; ++l)
                {
                    int prodIndex = i * mainElement.numNodes * uNum + j * uNum + l;

                    int vecIndex = i * mainElement.numNodes * uNum * 3 + k * uNum * 3 + l * 3;

                    prod[prodIndex] += mainElement.stiffnessMatrixX[stiffIndex] * flux.node[vecIndex] + \
                                       mainElement.stiffnessMatrixY[stiffIndex] * flux.node[vecIndex + 1] + \
                                       mainElement.stiffnessMatrixZ[stiffIndex] * flux.node[vecIndex + 2];

                }
            }

}