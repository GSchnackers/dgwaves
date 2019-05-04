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
                 const Properties & matProp, Quantity & flux){

    std::size_t i, j;

    // At the nodes.
    for(i = 0; i < u.node.size(); ++i)
    {
        int fluxIndex = 3 * i;
        int propIndex = i/6;

        switch (i % 6)
        {
            case 0: // Dx flux
                flux.node[fluxIndex] = 0;
                flux.node[fluxIndex + 1] = u.node[i + 5];
                flux.node[fluxIndex + 2] = -u.node[i + 4];
                break;

            case 1: // Dy flux
                flux.node[fluxIndex] = -u.node[i + 4];
                flux.node[fluxIndex + 1] = 0;
                flux.node[fluxIndex + 2] = u.node[i + 2];
                break;

            case 2: // Dz flux
                flux.node[fluxIndex] = u.node[i + 2];
                flux.node[fluxIndex + 1] = -u.node[i + 1];
                flux.node[fluxIndex + 2] = 0;
                break;

            case 3: // Bx flux
                flux.node[fluxIndex] = 0;
                flux.node[fluxIndex + 1] = u.node[i - 1];
                flux.node[fluxIndex + 2] = -u.node[i - 2];
                break;

            case 4: // By flux
                flux.node[fluxIndex] = -u.node[i - 2];
                flux.node[fluxIndex + 1] = 0;
                flux.node[fluxIndex + 2] = u.node[i - 4];
                break;

            case 5: // Bz flux
                flux.node[fluxIndex] = u.node[i - 4];
                flux.node[fluxIndex + 1] = -u.node[i - 5];
                flux.node[fluxIndex + 2] = 0;
                break;
        
        }

        if (!(i % 6) || i % 6 == 1 || i % 6 == 2)
            for(j = 0; j < 3; ++j)
                flux.node[fluxIndex + j] *= 1/matProp.relPermeability.node[propIndex];
        else
            for(j = 0; j < 3; ++j)
                flux.node[fluxIndex + j] *= 1/matProp.relPermittivity.node[propIndex];

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
                flux.gp[fluxIndex + j].first *= 1/matProp.relPermeability.gp[propIndex].first;
                flux.gp[fluxIndex + j].second *= 1/matProp.relPermeability.gp[propIndex].second;
            }
        else
            for(j = 0; j < 3; ++j)
            {
                flux.gp[fluxIndex + j].first *= 1/matProp.relPermittivity.gp[propIndex].first;
                flux.gp[fluxIndex + j].second *= 1/matProp.relPermittivity.gp[propIndex].second;
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

void numFluxELM(const Element & frontierElement, const Properties & matProp, const double alpha, \
                Quantity & flux){

    std::size_t i, j = 0;

    double factor1 = alpha * VACUUM_IMPEDANCE;
    double factor2 = alpha * VACUUM_CONDUCTANCE;

    if(alpha > 1 || alpha < 0)
    {
        gmsh::logger::write("The alpha coefficient is not suited for the Lax-Friedrichs method.", "error");
        exit(-1);
    }

    for(i = 0; i < flux.num.size(); ++i)
    {

        int propIndex = i/18;
        int bIndex = i + 3;

        if(!(i % 18) && i) j += 3;

        // Numerical flux associated to the electric fields.
        flux.num[i] = 1/(matProp.conductance.gp[propIndex].first + matProp.conductance.gp[propIndex].second) * \
                      (matProp.conductance.gp[propIndex].first * flux.gp[i].first + \
                      matProp.conductance.gp[propIndex].second * flux.gp[i].second \
                      + factor2 * (flux.gp[bIndex].first - flux.gp[bIndex].second)\
                      * frontierElement.normals[j + (bIndex % 3)]);

        // Numerical fluxes associated to magnetic fields.
        flux.num[bIndex] = 1/(matProp.impedance.gp[propIndex].first + matProp.impedance.gp[propIndex].second) * \
                           (matProp.impedance.gp[propIndex].first * flux.gp[bIndex].first + \
                           matProp.impedance.gp[propIndex].second * flux.gp[bIndex].second \
                           + factor1 * (flux.gp[i].first - flux.gp[i].second) \
                           * frontierElement.normals[j + (i % 3)]);

        if(i % 3 == 2) i += 3;

    }

}

