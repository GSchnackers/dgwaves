// This file contains all necessary functions to deal with the fluxes.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Simple flux where vec(f) = vec(c)u, where c is a constant vector representing a velocity. 
void physFluxCu(const Quantity & u, const Element & mainElement, const Element & frontierElement,\
                Quantity & flux){

    std::size_t i, j;
    std::vector<double> c = {1, 0, 0};

    flux.node.resize(mainElement.nodeTags.size() * 3);

    for(i = 0; i < mainElement.nodeTags.size(); ++i) // loop over the nodes of the main elements.
        for(j = 0; j < 3; ++j) // Loop over the components of the physical flux.
        {
            int index = i * 3 + j;
            flux.node[index] = c[j] * u.node[i];
        }

}

// Computes a fairly simple upwind flux: it choses the flux from which information comes from.
void numFluxCuUpwind(const Element & mainElement, Quantity & flux){

    std::size_t i, j;
    std::vector<double> c = {1, 0, 0};

    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
        {
            int scalarProd = 0;

            for(k = 0; k < 3; ++k)
            {
                int index = i * mainElement.numNodes * 3 + j * 3 + k;

                scalarProd += flux.node[index] * mainElement.normals[index];

            }

            if(scalarProd > 0)
                for(k = 0; k < 3; ++k)
                {
                    int index = i * frontierElement.numGp * 3 + j * 3 + k;

                    flux.numGp[index].first = flux.gp[index].first;
                    flux.numGp[index].second = -flux.numGp[index].first;

                }

            else if(scalarProd < 0)
                for(k = 0; k < 3; ++k)
                {
                    int index = i * frontierElement.numGp * 3 + j * 3 + k;

                    flux.numGp[index].first = flux.gp[index].second;
                    flux.numGp[index].second = -flux.numGp[index].first;

                }

            
        }    

}

