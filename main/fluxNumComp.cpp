// This functions computes the numerical fluxes at the gauss points at the frontier elements.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void fluxNumBackward(const Element & mainElement, const Element & frontierElement, \
                 const std::vector<double> & fluxPhysGp, const std::vector<double> & velocityGpFrontier,\
                 std::vector<double> & fluxNumGp){

    std::size_t i, j, k, l;

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numGp; ++j)
        {

            double scalarProd = 0;

            // Gets the indices of the neighbours of the frontier element.
            int neighbourIndex1 = frontierElement.neighbours[i].first;
            int neighbourIndex2 = frontierElement.neighbours[i].second;
            int fluxPhysRadical = neighbourIndex1;

            // The sign will play an important role 
            int sign = 1;

            for(k = 0; k < 3 ; ++k)
            {
                int index = i * frontierElement.numGp * 3 + j * 3 + k;
                scalarProd += velocityGpFrontier[index] * frontierElement.normals[index];
            }

            if(scalarProd < 0)
            {
                sign = -1;
                fluxPhysRadical = neighbourIndex2;
            }

            else if (scalarProd == 0) sign = 0;
            

            for(k = 0; k < mainElement.numGp; ++k)
                for(l = 0; l < 3; ++l)
                {
                    int fluxNumIndex1 = neighbourIndex1 * mainElement.numGp * 3 + k * 3 + l;
                    int fluxNumIndex2 = neighbourIndex2 * mainElement.numGp * 3 + k * 3 + l;
                    int fluxPhysIndex = fluxPhysRadical * mainElement.numGp * 3 + j * 3 + l;
                    
                    fluxNumGp[fluxNumIndex1] = sign * fluxPhysGp[fluxPhysIndex];
                    fluxNumGp[fluxNumIndex2] = - fluxNumGp[fluxNumIndex1];

                }

        }

}