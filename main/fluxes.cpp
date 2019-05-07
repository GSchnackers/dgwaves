// This file contains all necessary functions to deal with the fluxes.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Simple flux where vec(f) = vec(c)u, where c is a constant vector representing a velocity. 
void physFluxCu(const Quantity & u, const Element & mainElement, const Element & frontierElement,\
                Quantity & flux, std::vector<double> c){

    std::size_t i;

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
                 const Properties & matProp, Quantity & flux){

    std::size_t i, j;

    // At the nodes.
    for(i = 0; i < u.node.size(); ++i)
    {
        int fluxIndex = 3 * i;
        int propIndex = i/6;

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

        if (!(i % 6) || i % 6 == 1 || i % 6 == 2)
            for(j = 0; j < 3; ++j)
                flux.node[fluxIndex + j] *= -1/matProp.relPermittivity.node[propIndex];
        else
            for(j = 0; j < 3; ++j)
                flux.node[fluxIndex + j] *= 1/matProp.relPermeability.node[propIndex];

    }

    // At the gauss points.
    for(i = 0; i < u.gp.size(); ++i)
    {
        int fluxIndex = 3 * i;
        int propIndex = i/6;

        switch (i % 6)
        {
            case 0: // Ex flux

                flux.gp[fluxIndex].first = flux.gp[fluxIndex].second = 0;

                flux.gp[fluxIndex + 1].first = u.gp[i + 5].first;
                flux.gp[fluxIndex + 2].first = -u.gp[i + 4].first;
                
                flux.gp[fluxIndex + 1].second = u.gp[i + 5].second;
                flux.gp[fluxIndex + 2].second = -u.gp[i + 4].second;

                break;

            case 1: // Ey flux

                flux.gp[fluxIndex + 1].first = flux.gp[fluxIndex + 1].second = 0;;

                flux.gp[fluxIndex].first = -u.gp[i + 4].first;
                flux.gp[fluxIndex + 2].first = u.gp[i + 2].first;

                flux.gp[fluxIndex].second = -u.gp[i + 4].second;
                flux.gp[fluxIndex + 2].second = u.gp[i + 2].second;

                break;

            case 2: // Ez flux

                flux.gp[fluxIndex + 2].first = flux.gp[fluxIndex + 2].second = 0;

                flux.gp[fluxIndex].first = u.gp[i + 2].first;
                flux.gp[fluxIndex + 1].first = -u.gp[i + 1].first;

                flux.gp[fluxIndex].second = u.gp[i + 2].second;
                flux.gp[fluxIndex + 1].second = -u.gp[i + 1].second;
                

                break;

            case 3: // Hx flux

                flux.gp[fluxIndex].first = flux.gp[fluxIndex].second = 0;

                flux.gp[fluxIndex + 1].first = u.gp[i - 1].first;
                flux.gp[fluxIndex + 2].first = -u.gp[i - 2].first;

                flux.gp[fluxIndex + 1].second = u.gp[i - 1].second;
                flux.gp[fluxIndex + 2].second = -u.gp[i - 2].second;

                break;

            case 4: // Hy flux

                flux.gp[fluxIndex + 1].first = flux.gp[fluxIndex + 1].second = 0;

                flux.gp[fluxIndex].first = -u.gp[i - 2].first;
                flux.gp[fluxIndex + 2].first = u.gp[i - 4].first;

                flux.gp[fluxIndex].second = -u.gp[i - 2].second;
                flux.gp[fluxIndex + 2].second = u.gp[i - 4].second;

                break;

            case 5: // Hz flux

                flux.gp[fluxIndex + 2].first = flux.gp[fluxIndex + 2].second = 0;

                flux.gp[fluxIndex].first = u.gp[i - 4].first;
                flux.gp[fluxIndex + 1].first = -u.gp[i - 5].first;

                flux.gp[fluxIndex].second = u.gp[i - 4].second;
                flux.gp[fluxIndex + 1].second = -u.gp[i - 5].second;

                break;
        
        }

        if (!(i % 6) || i % 6 == 1 || i % 6 == 2)
            for(j = 0; j < 3; ++j)
            {
                flux.gp[fluxIndex + j].first *= -1/matProp.relPermittivity.gp[propIndex].first;
                flux.gp[fluxIndex + j].second *= -1/matProp.relPermittivity.gp[propIndex].second;
            }
        else
            for(j = 0; j < 3; ++j)
            {
                flux.gp[fluxIndex + j].first *= 1/matProp.relPermeability.gp[propIndex].first;
                flux.gp[fluxIndex + j].second *= 1/matProp.relPermeability.gp[propIndex].second;
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

// This function implements the lax-friedrichs flux for electromagnetic equations.
void numFluxELM(const Element & frontierElement, const double alpha, Quantity & u, Quantity & flux){

    std::size_t i;

    if(alpha > 1 || alpha < 0)
    {
        gmsh::logger::write("The alpha coefficient is not suited for the Lax-Friedrichs method.", "error");
        exit(-1);
    }

    for(i = 0; i < flux.num.size(); ++i)
        flux.num[i] = 0.5 * (flux.gp[i].first + flux.gp[i].second + \
                      alpha * (u.gp[i/3].first - u.gp[i/3].second) * frontierElement.normals[i/18 + (i % 3)]);


}

