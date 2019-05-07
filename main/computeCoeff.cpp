/*

 This file contains the function that computes the coefficients for the forward EUler method or for 
 Runge-Kutta 4. For more information, please see the documentation of each called function within.

*/

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void computeCoeff(const Element & mainElement, const Element & frontierElement,\
                  const std::vector<Parameter> & bcParam, const Simulation & simulation,\
                  const Properties & matProp, const double t, Quantity & u, Quantity & flux, \
                  std::vector<double> & k){

    std::vector<double> SFProd(simulation.uNum * mainElement.nodeTags.size(), 0);
    std::vector<double> fluxVector(simulation.uNum * mainElement.nodeTags.size(), 0);

    computeBoundaryCondition(u, t, bcParam);
    valGp(u, mainElement, frontierElement, simulation.uNum);

    if(simulation.uNum == 6)
    {
        physFluxELM(u, frontierElement, mainElement, matProp, flux);
        numFluxELM(frontierElement, simulation.alpha, u, flux);
    }

    else if(simulation.uNum == 1)
    {
        std::vector<double> c = {1 , 0 , 0};
        physFluxCu(u, mainElement, frontierElement, flux, c);
        numFluxUpwind(frontierElement, flux);
    }

    stiffnessFluxProd(mainElement, flux, SFProd, simulation.uNum);
    std::cout << "Hello" << std::endl;
    numFluxIntegration(flux, mainElement, frontierElement, fluxVector, simulation.uNum);
    timeMarching(mainElement, SFProd, fluxVector, k, simulation.uNum);
    std::cout << "hello" << std::endl;

    if(simulation.debug) timeChecker(mainElement, frontierElement, flux, u, SFProd, fluxVector, t, simulation.uNum);
        

}