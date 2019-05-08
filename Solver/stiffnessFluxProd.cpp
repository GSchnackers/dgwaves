#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "structures.hpp"

// This function computes the product of the stiffness matrices and the physical flux vector at nodal values.
void stiffnessFluxProd(const Element & mainElement, const Quantity & flux, std::vector<double> & prod, int uNum){

    std::size_t i, j, k, l;

    // Components along x, y and z of the gradient.
    std::vector<double> fx(flux.node.size()/3, 0), fy(flux.node.size()/3, 0), fz(flux.node.size()/3, 0);
    std::fill(prod.begin(), prod.end(), 0);

    for(i = 0; i < flux.node.size(); i += 3)
    {
        fx[i/3] = flux.node[i];
        fy[i/3] = flux.node[i + 1];
        fz[i/3] = flux.node[i + 2];
    }
    
    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < mainElement.numNodes; ++k)
            {
                int stiffIndex = i * mainElement.numNodes * mainElement.numNodes + j * mainElement.numNodes + k;

                for(l = 0; l < uNum; ++l)
                {
                    int prodIndex = i * mainElement.numNodes * uNum + j * uNum + l;

                    int vecIndex = i * mainElement.numNodes * uNum + k * uNum + l;

                    prod[prodIndex] += mainElement.stiffnessMatrixX[stiffIndex] * fx[vecIndex] + \
                                       mainElement.stiffnessMatrixY[stiffIndex] * fy[vecIndex] + \
                                       mainElement.stiffnessMatrixZ[stiffIndex] * fz[vecIndex];

                }
            }

}