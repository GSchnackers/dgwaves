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
                  const std::vector<Parameter> & bcParam, const double simStep, const double t, \
                  const Quantity & impedance, Quantity & u, Quantity & flux, std::vector<double> & k, \
                  const int debug, const double alpha){

    std::vector<double> SFProd(mainElement.nodeTags.size(), 0);
    std::vector<double> fluxVector(mainElement.nodeTags.size(), 0);

    computeBoundaryCondition(u, t, bcParam);
    valGp(u, mainElement, frontierElement);
    physFluxELM(u, frontierElement, mainElement, flux);
    numFluxELM(frontierElement, impedance, alpha, flux); 
    stiffnessFluxProd(mainElement, flux, SFProd);
    numFluxIntegration(flux, mainElement, frontierElement, fluxVector);
    timeMarching(mainElement, SFProd, fluxVector, k);

    if(debug) timeChecker(mainElement, frontierElement, flux, u, SFProd, fluxVector, t);
        

}