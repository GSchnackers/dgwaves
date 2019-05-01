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

    std::vector<double> SFProd(6 * mainElement.nodeTags.size(), 0);
    std::vector<double> fluxVector(6 * mainElement.nodeTags.size(), 0);

    computeBoundaryCondition(u, t, bcParam);
    valGp(u, mainElement, frontierElement, 6);
    physFluxELM(u, frontierElement, mainElement, matProp, flux);
    numFluxELM(frontierElement, matProp, simulation.alpha, flux);
    stiffnessFluxProd(mainElement, flux, SFProd);
    numFluxIntegration(flux, mainElement, frontierElement, fluxVector);
    timeMarching(mainElement, SFProd, fluxVector, k);

    if(simulation.debug) timeChecker(mainElement, frontierElement, flux, u, SFProd, fluxVector, t);
        

}