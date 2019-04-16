// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(Element & mainElement, Element & frontierElement, View & mainView){

    double t, step = 0.001; // t is the time.
    std::size_t i, j, k; // loop variables.

    Quantity u; // unknowns of the problem.
    Quantity flux; // fluxs of the problem.
    
    std::vector<double> SFProd;
    std::vector<double> fluxVector;

    // Initialization of the nodal values.
    u.node.resize(mainElement.nodeTags.size(), 0);
    u.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    u.bound.resize(mainElement.nodeTags.size(), 0);
    u.boundSign.resize(mainElement.nodeTags.size(), 0);

    flux.node.resize(mainElement.nodeTags.size() * 3, 0);
    flux.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp * 3, std::make_pair(0,0));
    flux.direction.resize(flux.gp.size(), 0);
    flux.num.resize(flux.gp.size(), 0);
    flux.bound.resize(mainElement.nodeTags.size() * 3, 0);

    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type...");
    setBoundaryConditions(mainElement, u);
    std::cout << "Done." << std::endl;

    SFProd.resize(mainElement.nodeTags.size(), 0);
    fluxVector.resize(mainElement.nodeTags.size(), 0);

    for(t = 0; t < 0.007; t += step)
    {    
        computeBoundaryCondition(mainElement, u, t);

        // Boundary Conditions verification.

        std::cout << "BC's verifier at t = " << t << std::endl;
        for(i = 0; i < u.bound.size(); ++i)
            std::cout << "Element: " << mainElement.elementTag[i/mainElement.numNodes] << " Node: " << mainElement.nodeTags[i] << " Value: " << u.bound[i] << " " << std::endl;
        std::cout << std::endl;
        
        valGp(u, mainElement, frontierElement, 1);
        std::cout << "Gauss points verifier at t = " << t << std::endl;
        
        for(i = 0; i < u.gp.size(); ++i)
        {
            std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/frontierElement.numGp].first] << " Gauss Point: " << i % frontierElement.numGp << " Value: " << u.gp[i].first;
            if(frontierElement.neighbours[i/frontierElement.numGp].second >= 0)
                std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/frontierElement.numGp].second] << " Gauss Point: " << i % frontierElement.numGp << " Value: " << u.gp[i].second << " ";
            else
                std::cout << " / Element: NONE Gauss Point: " << i % frontierElement.numGp << " Value: " << u.gp[i].second << " ";    
            std::cout << std::endl;
        }
        std::cout << std::endl;
        /*physFluxCu(u, mainElement, frontierElement, flux);
        
        numFluxUpwind(frontierElement, flux);
        stiffnessFluxProd(mainElement, flux, SFProd);
        numFluxIntegration(flux, mainElement, frontierElement, fluxVector);

       /*  for(i = 0; i < fluxVector.size(); ++i)
            std::cout << fluxVector[i] << " " << SFProd[i] << " " << mainElement.nodeTags[i] << std::endl;

        for(i = 0; i < mainElement.elementTag.size(); ++i)
        {
            std::cout << "Element " << mainElement.elementTag[i] << std::endl;
            for(j = 0; j < mainElement.numNodes; ++j)
                std::cout << fluxVector[i * mainElement.numNodes + j] << std::endl;

            std::cout << std::endl;

        } */
        
        /* for(i = 0; i < mainElement.elementTag.size(); ++i)
            for(j = 0; j < mainElement.numNodes; ++j)
            {
                int uIndex = i * mainElement.numNodes + j;
                double tmpProd = 0;

                for(k = 0; k < mainElement.numNodes; ++k)
                {
            
                    int matrixIndex = i * mainElement.numNodes * mainElement.numNodes + \
                                      j * mainElement.numNodes + k;

                    int vecIndex = i * mainElement.numNodes + k;

                    tmpProd += mainElement.massMatrixInverse[matrixIndex] * \
                                       (SFProd[vecIndex] - fluxVector[vecIndex]);
                    //std::cout << fluxVector[vecIndex] << std::endl;
                    
                }


                mainView.data[i][j] = u.node[uIndex] += step * tmpProd;

                //std::cout << u.node[uIndex] << " " << mainElement.nodeTags[uIndex] << std::endl;
                

            }*/   
            
        gmsh::view::addModelData(mainView.tag, int(t/step), mainView.modelName, mainView.dataType, \
                                 mainElement.elementTag, mainView.data, t, 1);
        

    }

    
    gmsh::view::write(mainView.tag, "results.msh");

}