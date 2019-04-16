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

    for(t = 0; t < 0.01; t += step)
    {    
        computeBoundaryCondition(mainElement, u, t);

        // Boundary Conditions verification.
        std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << " %%%%%%%%%%%%%%%%%%%%% TIME STEP t = " << t << " %%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

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

        physFluxCu(u, mainElement, frontierElement, flux);

        std::cout << "Physical nodal flux verifier at t = " << t << std::endl;
        for(i = 0; i < flux.node.size(); ++i)
            std::cout << "Element: " << mainElement.elementTag[i/(3 * mainElement.numNodes)] << " Node: " << mainElement.nodeTags[i/3] << " Value: " << flux.node[i] << std::endl;
        std::cout << std::endl;

        std::cout << "Physical flux at the Gauss points verifier at t = " << t << std::endl;
        
        for(i = 0; i < flux.gp.size(); ++i)
        {
            std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 3)].first] << " Gauss Point: " << (i / 3) % frontierElement.numGp << " Value: " << flux.gp[i].first;
            if(frontierElement.neighbours[i/(3*frontierElement.numGp)].second >= 0)
                std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 3)].second] << " Gauss Point: " << (i / 3) % frontierElement.numGp << " Value: " << flux.gp[i].second << " ";
            else
                std::cout << " / Element: NONE Gauss Point: " << (i / 3) % frontierElement.numGp << " Value: " << flux.gp[i].second << " ";    
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "Physical flux direction at the Gauss points at t = " << t << std::endl;
        for(i = 0; i < flux.direction.size(); ++i)
        {
            std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 3)].first];
            if(frontierElement.neighbours[i/(3*frontierElement.numGp)].second >= 0)
                std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 3)].second] << " Gauss Point: " << (i / 3) % frontierElement.numGp << " Value: " << flux.direction[i] << " ";
            else
                std::cout << " / Element: NONE Gauss Point: " << (i / 3) % frontierElement.numGp << " Value: " << flux.direction[i] << " ";    
            std::cout << std::endl;
        }
        std::cout << std::endl;
        
        numFluxUpwind(frontierElement, flux);
        std::cout << "Numerical flux at the Gauss points verifier at t = " << t << std::endl;
        for(i = 0; i < flux.num.size(); ++i)
        {
            std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 3)].first];
            if(frontierElement.neighbours[i/(3*frontierElement.numGp)].second >= 0)
                std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 3)].second] << " Gauss Point: " << (i / 3) % frontierElement.numGp << " Value: " << flux.num[i] << " ";
            else
                std::cout << " / Element: NONE Gauss Point: " << (i / 3) % frontierElement.numGp << " Value: " << flux.num[i] << " ";    
            std::cout << std::endl;
        }
        std::cout << std::endl;        

        stiffnessFluxProd(mainElement, flux, SFProd);
        std::cout << "Stiffness/Flux product verifier at t = " << t << std::endl;
        for(i = 0; i < SFProd.size(); ++i)
            std::cout << "Element: " << mainElement.elementTag[i/mainElement.numNodes] << " Node: " << mainElement.nodeTags[i] << " Value: " << SFProd[i] << std::endl;
        std::cout << std::endl;
        /*numFluxIntegration(flux, mainElement, frontierElement, fluxVector);

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