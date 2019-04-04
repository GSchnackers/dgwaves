// This functions computes the physical flux at the gauss points of the elements in parametric coordinates and
// the physical flux of the nodes in real coordinates.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Remark: here we simply use as a first test a physical flux of the form c u, where c is the velocity of
// the propagated quantity u.

void fluxComp(const Element & element, const std::vector<double> u, const std::vector<double> & uGp,\
              const std::vector<double> & velocity, const std::vector<double> & velocityGp,\
              std::vector<double> & fluxPhysGpParam, std::vector<double> & fluxPhys){

    std::size_t i, j, k, l;

    for(i = 0; i < element.elementTag.size(); ++i)
        for(j = 0; j < 3; ++j)
            for(k = 0; k < element.numNodes; ++k)
            {
                int nodeIndex = i * element.numNodes * 3 + k * 3 + j;

                fluxPhys[nodeIndex] = u[nodeIndex/3] * velocity[nodeIndex];

                for(l = 0; l < element.numGp; ++l)
                {
                    int gpIndex = i * element.numGp * 3 + l * 3 + j;
                    fluxPhysGpParam[gpIndex] = uGp[gpIndex/3] + velocityGp[gpIndex];
                }

            }

}