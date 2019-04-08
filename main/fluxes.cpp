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

    flux.node.resize(mainElement.nodeTags.size() * 3);
    flux.numGp.resize(frontierElement.elementTag.size() * frontierElement.numGp * 3);

    for(i = 0; i < mainElement.nodeTags.size(); ++i) // loop over the nodes of the main elements.
        for(j = 0; j < 3; ++j) // Loop over the components of the physical flux.
        {
            int index = i * 3 + j;
            flux.node[index] = c[j] * u.node[i];
        }


    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numGp; ++j)
            for(k = 0; k < 3; ++k)
            {
                int index = i * frontierElement.numGp * 3 + j * 3 + k;
                flux.numGp[index].first = u.numGp[index/3].first * c[k];
                flux.numGp[index].second = u.numGp[index/3].second * c[k];
            }

}

// Computes a fairly simple upwind flux: it choses the flux from which information comes from.
void numFluxUpwind(const Element & frontierElement, Quantity & flux){

    std::size_t i, j, k;

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
    {
        for(j = 0; j < frontierElement.numGp; ++j)
        {
            int scalarProd = 0;
            int fluxIndex = i * frontierElement.numGp * 3 + j * 3;

            // Computation of the scalar product at the gauss point.
            for(k = 0; k < 3; ++k)
            {
                int index = fluxIndex + k;
                scalarProd += flux.numGp[index].first * frontierElement.normals[index];
            }

            // Selection of the upwind flux.
            if(scalarProd > 0)
                for(k = 0; k < 3; ++k)
                {
                    int index = fluxIndex + k;
                    flux.numGp[index].second = -flux.numGp[index].first;
                }

            else if(scalarProd <= 0)   
                 for(k = 0; k < 3; ++k)
                 {
                    int index = fluxIndex + k;
                    flux.numGp[index].first = -flux.numGp[index].second;
                 }
            
            
        } 
    }   

}

