// This file contains all necessary functions to deal with the fluxes.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Simple flux where vec(f) = vec(c)u, where c is a constant vector representing a velocity. 
void physFluxCu(const Quantity & u, const Element & mainElement, const Element & frontierElement,\
                Quantity & flux){

    std::size_t i, j, k;
    std::vector<double> c = {1, 0, 0};

    for(i = 0; i < mainElement.nodeTags.size(); ++i) // loop over the nodes of the main elements.
        for(j = 0; j < 3; ++j) // Loop over the components of the physical flux.
        {
            flux.node[i * 3 + j] = c[j] * u.node[i];
            std::cout << flux.node[i * 3 + j] << std::endl;
        }
        

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numGp; ++j)
        {
            int smallIndex = i * frontierElement.numGp + j;
            
            for(k = 0; k < 3; ++k)
            {
                int index = i * frontierElement.numGp * 3 + j * 3 + k;

                flux.gp[index].first = u.gp[smallIndex].first * c[k];
                flux.gp[index].second = u.gp[smallIndex].second * c[k];


                flux.direction[index] = c[k];
            }
        }

}

// Computes a fairly simple upwind flux: it choses the flux from which information comes from.
void numFluxUpwind(const Element & frontierElement, Quantity & flux){

    std::size_t i, j, k;

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
    {
        for(j = 0; j < frontierElement.numGp; ++j)
        {
            double scalarProd = 0;
            int fluxIndex = i * frontierElement.numGp * 3 + j * 3;

            // Computation of the scalar product at the gauss point.
            for(k = 0; k < 3; ++k)
            {
                int index = fluxIndex + k;
                scalarProd += flux.direction[index] * frontierElement.normals[index];
            }

            // Selection of the upwind flux.
            
            for(k = 0; k < 3; ++k)
            {
                int index = fluxIndex + k;

                if(scalarProd > 0)
                    flux.num[index] = flux.gp[index].first;

                else if(scalarProd < 0)
                    flux.num[index] = flux.gp[index].second;
                else
                    flux.num[index] = 0;

            }
            
        } 

    } 

}

