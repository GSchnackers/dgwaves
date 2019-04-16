// This file contains all necessary functions to deal with the fluxes.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Simple flux where vec(f) = vec(c)u, where c is a constant vector representing a velocity. 
void physFluxCu(const Quantity & u, const Element & mainElement, const Element & frontierElement,\
                Quantity & flux){

    std::size_t i;
    std::vector<double> c = {1 , 0 , 0};

    for(i = 0; i < flux.node.size(); ++i) // loop over the nodes of the main elements.
        flux.node[i] = c[i % 3] * u.node[i/3];
    
    
    for(i = 0; i < flux.gp.size(); ++i) 
    {
        flux.gp[i].first = c[i % 3] * u.gp[i/3].first;
        flux.gp[i].second = c[i % 3] * u.gp[i/3].second;
        flux.direction[i] = c[i % 3];
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

            // Computation of the scalar product at the gauss point.
            for(k = 0; k < 3; ++k)
            {
                int index = i * frontierElement.numGp * 3 + j * 3 + k;
                scalarProd += flux.direction[index] * frontierElement.normals[index];
            }

            // Selection of the upwind flux.
            
            for(k = 0; k < 3; ++k)
            {
                int index = i * frontierElement.numGp * 3 + j * 3 + k;

                if(scalarProd > 0)
                    flux.num[index] = flux.gp[index].first;

                else if(scalarProd = 0)
                    flux.num[index] = flux.gp[index].second;
                    
                else
                    flux.num[index] = 0;

            }
            
        } 

    } 

}

