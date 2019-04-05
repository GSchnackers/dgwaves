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
    flux.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp * 3);

    for(i = 0; i < mainElement.nodeTags.size(); ++i) // loop over the nodes of the main elements.
        for(j = 0; j < 3; ++j) // Loop over the components of the physical flux.
        {
            int index = i * 3 + j;
            flux.node[index] = c[j] * u.node[i];
        }
    
    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numNodes; ++j)
            for(k = 0; k < 3; ++k)
            {
                int index = i * frontierElement.numNodes * 3 + j * 3 + k;
                
                flux.gp[index].first = u.gp[index].first * c[k];
                flux.gp[index].second = u.gp[index].second * c[k];
            }

}

