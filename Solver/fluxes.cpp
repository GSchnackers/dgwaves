// This file contains all necessary functions to deal with the fluxes.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "structures.hpp"
#include <omp.h>

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
                 const Properties & matProp, double t, Quantity & flux){

    std::size_t i, j;

    // At the nodes.
    #pragma omp parallel for default(shared) private(j)
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
    #pragma omp parallel for default(shared) private(j)
    for(i = 0; i < u.gp.size(); ++i)
    {
        int fluxIndex = 3 * i;
        int propIndex = i/6;
        bool condition = frontierElement.neighbours[i/(6 * frontierElement.numGp)].second >= 0;
    
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
                
                flux.gp[fluxIndex + 1].first = flux.gp[fluxIndex + 1].second = 0;
                
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

        if ((i % 6) < 3)
            for(j = 0; j < 3; ++j)
            {
                flux.gp[fluxIndex + j].first *= -1/matProp.relPermittivity.gp[propIndex].first;

                if(condition)
                    flux.gp[fluxIndex + j].second *= -1/matProp.relPermittivity.gp[propIndex].second;
                else
                    flux.gp[fluxIndex + j].second *= -1/matProp.relPermittivity.gp[propIndex].first;
            }
            
        else
            for(j = 0; j < 3; ++j)
            {
                flux.gp[fluxIndex + j].first *= 1/matProp.relPermeability.gp[propIndex].first;

                if(condition)
                    flux.gp[fluxIndex + j].second *= 1/matProp.relPermeability.gp[propIndex].second;
                else
                    flux.gp[fluxIndex + j].second *= 1/matProp.relPermeability.gp[propIndex].first;
                
            }
        

    }

}

// Computes a fairly simple upwind flux: it choses the flux from which information comes from.
void numFluxUpwind(const Element & frontierElement, Quantity & flux){

    std::size_t i, j, k;

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
    {
        
    }
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
                    flux.num[index/3] += flux.gp[index].first * frontierElement.normals[index];

                else if(scalarProd < 0)
                    flux.num[index/3] += flux.gp[index].first * frontierElement.normals[index];

            }
            
        } 

    } 

}

// This function implements the numerical flux for electromagnetic equations.
void numFluxELM(const Element & frontierElement, const Properties & matProp, const double alpha, Quantity & u, \
                Quantity & flux){

    std::size_t i, j;
    std::vector<double> physFluxNorm(flux.gp.size()/3, 0);
    double scalU; 

    if(alpha > 1 || alpha < 0)
    {
        gmsh::logger::write("The alpha coefficient is not suited for the Lax-Friedrichs method.", "error");
        exit(-1);
    }

    // scalar product between the normal and the flux gauss point.
    for(i = 0; i < flux.gp.size(); ++i)
        physFluxNorm[i/3] += frontierElement.normals[i/18 * 3 + (i % 3)] * \
                             (matProp.speedGp[i/18].first * flux.gp[i].first + matProp.speedGp[i/18].second * flux.gp[i].second);

    #pragma omp parallel for default(shared) private(j, scalU)
    for(i = 0; i < flux.num.size(); ++i)
    {
        if(!(i % 6) || i % 6 == 3)
        {
            scalU = 0;

            for(j = 0; j < 3; ++j)
                scalU += (u.gp[i + j].first - u.gp[i + j].second) * frontierElement.normals[i/6 * 3 + j];

        }
        
        flux.num[i] = matProp.speedGpSumInv[i/6] * (physFluxNorm[i] + alpha * (u.gp[i].first - u.gp[i].second) - scalU * frontierElement.normals[i/6 * 3 + (i % 3)]);
    }


}