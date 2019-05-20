/*
    This file deals with the computation of the coefficients for the 
*/

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "solver.hpp"
#include "structures.hpp"

// Function that check all relevant quantities over time if required. It is used only in debug mode.
void timeChecker(const Element & mainElement, const Element & frontierElement,\
                 const Quantity & flux, const Quantity & u, const std::vector<double> & SFProd, \
                 const std::vector<double> & fluxVector, const double t, int numU){

    std::size_t i, j, k;

    // Boundary Conditions verification.
    std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << " %%%%%%%%%%%%%%%%%%%%% TIME STEP t = " << t << " %%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    std::cout << "Gauss points verifier at t = " << t << std::endl;
        
    // Value of u at the gauss Point verifier.
    for(i = 0; i < u.gp.size(); ++i)
    {
        std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * numU)].first] << " Gauss Point: " << i % (frontierElement.numGp * numU) << " Value: " << u.gp[i].first;
        if(frontierElement.neighbours[i/(frontierElement.numGp * numU)].second >= 0)
            std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * numU)].second] << " Gauss Point: " << i % (frontierElement.numGp * numU) << " Value: " << u.gp[i].second << " ";
        else
            std::cout << " / Element: NONE Gauss Point: " << i % (frontierElement.numGp*numU) << " Value: " << u.gp[i].second << " ";    
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Value of the nodal flux verifier.
    std::cout << "Physical nodal flux verifier at t = " << t << std::endl;
    for(i = 0; i < flux.node.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(3 * numU * mainElement.numNodes)] << " Node: " << mainElement.nodeTags[i/(3 * numU)] << " Value: " << flux.node[i] << std::endl;
    std::cout << std::endl;

    std::cout << "Physical flux at the Gauss points verifier at t = " << t << std::endl;
    
    // Value of the flux at the gauss point verifyer.
    for(i = 0; i < flux.gp.size(); ++i)
    {
        std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 3 * numU)].first] << " Gauss Point: " << (i / (3 * numU)) % frontierElement.numGp << " Value: " << flux.gp[i].first;
        if(frontierElement.neighbours[i/((3 * numU)*frontierElement.numGp)].second >= 0)
            std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * (3 * numU))].second] << " Gauss Point: " << (i / (3 * numU)) % frontierElement.numGp << " Value: " << flux.gp[i].second << " ";
        else
            std::cout << " / Element: NONE Gauss Point: " << (i / (3 * numU)) % frontierElement.numGp << " Value: " << flux.gp[i].second << " ";    
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Value of the numerical flux verifier.
    std::cout << "Numerical flux at the Gauss points verifier at t = " << t << std::endl;
    for(i = 0; i < flux.num.size(); ++i)
    {
        std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * (3 * numU))].first];
        if(frontierElement.neighbours[i/((3 * numU)*frontierElement.numGp)].second >= 0)
            std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * (3 * numU))].second] << " Gauss Point: " << (i / (3 * numU)) % frontierElement.numGp << " Value: " << flux.num[i] << " ";
        else
            std::cout << " / Element: NONE Gauss Point: " << (i / (3 * numU)) % frontierElement.numGp << " Value: " << flux.num[i] << " ";    
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Stiffness flux product verifier.
    std::cout << "Stiffness/Flux product verifier at t = " << t << std::endl;
    for(i = 0; i < SFProd.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(mainElement.numNodes * numU)] << " Node: " << mainElement.nodeTags[i/numU] << " Value: " << SFProd[i] << std::endl;
    std::cout << std::endl;

    // Flux Integration verifier.
    std::cout << "Flux integration verifier at t = " << t << std::endl;
    for(i = 0; i < SFProd.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(mainElement.numNodes * numU)] << " Node: " << mainElement.nodeTags[i/numU] << " Value: " << fluxVector[i] << std::endl;
    std::cout << std::endl;

    // Nodal values verifier.
    std::cout << "Nodal value verifier at t = " << t << std::endl;
    for(i = 0; i < u.node.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(mainElement.numNodes * numU)] << " Node: " << mainElement.nodeTags[i/numU] << " Value: " << u.node[i] << std::endl;
    std::cout << std::endl;

}

void timeMarching(const Element & mainElement, const std::vector<double> & SFProd, \
                  const std::vector<double> & fluxVector, std::vector<double> & kVector, int uNum){

    std::size_t i, j, k, l;
    std::fill(kVector.begin(), kVector.end(), 0);

    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < mainElement.numNodes; ++k)
                for(l = 0; l < uNum; ++l)
                {
                    int uIndex      = i * mainElement.numNodes * uNum + j * uNum + l;
                    int matrixIndex = i * mainElement.numNodes * mainElement.numNodes + \
                                      j * mainElement.numNodes + k;

                    int vecIndex = i * mainElement.numNodes * uNum + k * uNum + l;

                    kVector[uIndex] += mainElement.massMatrixInverse[matrixIndex] * \
                                       (SFProd[vecIndex] - fluxVector[vecIndex]);
                    
                }

}

void computeCoeff(const Element & mainElement, const Element & frontierElement, const Simulation & simulation,\
                  const Properties & matProp, const double t, Quantity & u, Quantity & flux, \
                  std::vector<double> & k){

    std::vector<double> SFProd(simulation.uNum * mainElement.nodeTags.size(), 0);
    std::vector<double> fluxVector(simulation.uNum * mainElement.nodeTags.size(), 0);

    valGp(u, mainElement, frontierElement, simulation.uNum, matProp, t);

    if(simulation.uNum == 6)
    {
        physFluxELM(u, frontierElement, mainElement, matProp, t, flux);
        numFluxELM(frontierElement, matProp, simulation.alpha, u, flux);
    }

    else if(simulation.uNum == 1)
    {
        physFluxCu(u, mainElement, frontierElement, flux, simulation.c);
        numFluxUpwind(frontierElement, u, flux, simulation.alpha);
    }

    stiffnessFluxProd(mainElement, flux, SFProd, simulation.uNum);
    numFluxIntegration(flux, mainElement, frontierElement, fluxVector, simulation.uNum);
    timeMarching(mainElement, SFProd, fluxVector, k, simulation.uNum);

    if(simulation.debug) timeChecker(mainElement, frontierElement, flux, u, SFProd, fluxVector, t, simulation.uNum);
        

}