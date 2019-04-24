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

void physFluxELM(const Quantity & u, const Element & frontierElement, const Element & mainElement,\
                 Quantity & flux){

    std::size_t i, j, k;

    // At the nodes.
    for(i = 0; i < u.node.size(); ++i)
    {
        int fluxIndex = 3 * i;

        switch (i % 6)
        {
            case 0: // Ex flux
                flux.node[fluxIndex] = 0;
                flux.node[fluxIndex + 1] = u.node[i + 5];
                flux.node[fluxIndex + 2] = -u.node[i + 4];
                break;

            case 1: // Ey flux
                flux.node[fluxIndex] = -u.node[i + 4];
                flux.node[fluxIndex + 1] = 0;
                flux.node[fluxIndex + 2] = u.node[i + 2];
                break;

            case 2: // Ez flux
                flux.node[fluxIndex] = u.node[i + 2];
                flux.node[fluxIndex + 1] = -u.node[i + 1];
                flux.node[fluxIndex + 2] = 0;
                break;

            case 3: // Hx flux
                flux.node[fluxIndex] = 0;
                flux.node[fluxIndex + 1] = u.node[i - 1];
                flux.node[fluxIndex + 2] = -u.node[i - 2];
                break;

            case 4: // Hy flux
                flux.node[fluxIndex] = -u.node[i - 2];
                flux.node[fluxIndex + 1] = 0;
                flux.node[fluxIndex + 2] = u.node[i - 4];
                break;

            case 5: // Hz flux
                flux.node[fluxIndex] = u.node[i - 4];
                flux.node[fluxIndex + 1] = -u.node[i - 5];
                flux.node[fluxIndex + 2] = 0;
                break;
        
        }

    }


    // At the gauss points.
    for(i = 0; i < u.gp.size(); ++i)
    {
        int fluxIndex = 3 * i;
        int neighIndex1 = 6 * (frontierElement.neighbours[i].first * mainElement.numNodes + \
                          frontierElement.nodeCorrespondance[i].first);
        int neighIndex2 = 6 * (frontierElement.neighbours[i].second * mainElement.numNodes + \
                          frontierElement.nodeCorrespondance[i].second);

        switch (i % 6)
        {
            case 0: // Ex flux

                flux.gp[fluxIndex].first = 0;
                flux.gp[fluxIndex + 1].first = u.node[neighIndex1 + 5];
                flux.gp[fluxIndex + 2].first = -u.node[neighIndex1 + 4];

                flux.gp[fluxIndex].second = 0;
                flux.gp[fluxIndex + 1].second = u.node[neighIndex2 + 5];
                flux.gp[fluxIndex + 2].second = -u.node[neighIndex2 + 4];

                break;

            case 1: // Ey flux

                flux.gp[fluxIndex].first = -u.node[neighIndex1 + 4];
                flux.gp[fluxIndex + 1].first = 0;
                flux.gp[fluxIndex + 2].first = u.node[neighIndex1 + 2];

                flux.gp[fluxIndex].first = -u.node[neighIndex2 + 4];
                flux.gp[fluxIndex + 1].first = 0;
                flux.gp[fluxIndex + 2].first = u.node[neighIndex2 + 2];

                break;

            case 2: // Ez flux

                flux.gp[fluxIndex].first = u.node[neighIndex1 + 2];
                flux.gp[fluxIndex + 1].first = -u.node[neighIndex1 + 1];
                flux.gp[fluxIndex + 2].first = 0;

                flux.gp[fluxIndex].second = u.node[neighIndex2 + 2];
                flux.gp[fluxIndex + 1].second = -u.node[neighIndex2 + 1];
                flux.gp[fluxIndex + 2].second = 0;

                break;

            case 3: // Hx flux

                flux.gp[fluxIndex].first = 0;
                flux.gp[fluxIndex + 1].first = u.node[neighIndex1 - 1];
                flux.gp[fluxIndex + 2].first = -u.node[neighIndex1 - 2];

                flux.gp[fluxIndex].second = 0;
                flux.gp[fluxIndex + 1].second = u.node[neighIndex2 - 1];
                flux.gp[fluxIndex + 2].second = -u.node[neighIndex2 - 2];

                break;

            case 4: // Hy flux

                flux.gp[fluxIndex].first = -u.node[neighIndex1 - 2];
                flux.gp[fluxIndex + 1].first = 0;
                flux.gp[fluxIndex + 2].first = u.node[neighIndex1 - 4];

                flux.gp[fluxIndex].second = -u.node[neighIndex2 - 2];
                flux.gp[fluxIndex + 1].second = 0;
                flux.gp[fluxIndex + 2].second = u.node[neighIndex2 - 4];

                break;

            case 5: // Hz flux

                flux.gp[fluxIndex].first = u.node[neighIndex1 - 4];
                flux.gp[fluxIndex + 1].first = -u.node[neighIndex1 - 5];
                flux.gp[fluxIndex + 2].first = 0;

                flux.gp[fluxIndex].second = u.node[neighIndex2 - 4];
                flux.gp[fluxIndex + 1].second = -u.node[neighIndex2 - 5];
                flux.gp[fluxIndex + 2].second = 0;

                break;
        
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

                else if(scalarProd < 0)
                    flux.num[index] = flux.gp[index].second;
                    
                else
                    flux.num[index] = 0;

            }
            
        } 

    } 

}

void numFluxELM(const Element & frontierElement, const Quantity & impedance, const double alpha, Quantity & flux){

    std::size_t i, j, k;

    if(alpha > 1 || alpha < 0)
    {
        gmsh::logger::write("The alpha coefficient is not suited for the Lax-Friedrichs method.", "error");
        exit(-1);
    }

    for(i = 0; i < flux.num.size(); ++i)
    {

        int impIndex = i/18;

        // Numerical flux associated to the electric fields.
        flux.num[i] = (impedance.gp[impIndex].first + impedance.gp[impIndex].second)/2 * \
                      (1/impedance.gp[impIndex].first * flux.gp[i].first + \
                       1/impedance.gp[impIndex].second * flux.gp[i].second)\
                      + alpha/2 * (flux.gp[i + 3].first - flux.gp[i + 3].second);

        // Numerical fluxes associated to magnetic fields.
        flux.num[i + 3] = 2/(impedance.gp[impIndex + 3].first + impedance.gp[impIndex + 3].second) * \
                      (impedance.gp[impIndex + 3].first * flux.gp[i + 3].first + \
                       impedance.gp[impIndex + 3].second * flux.gp[i + 3].second)\
                      + alpha/2 * (flux.gp[i].first - flux.gp[i].second) * frontierElement.normals[i];

        if(i % 3 == 2) i += 3;

    }

}

