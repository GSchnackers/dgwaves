/*

 This file contains the function that computes the coefficients for the forward EUler method or for 
 Runge-Kutta 4. For more information, please see the documentation of each called function within.

*/

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void computeCoeff(const Element & mainElement, const Element & frontierElement, const double step, \
                  const double t, Quantity & u, Quantity & flux, std::vector<double> & k, int debug){

    std::vector<double> SFProd(mainElement.nodeTags.size(), 0);
    std::vector<double> fluxVector(mainElement.nodeTags.size(), 0);

    computeBoundaryCondition(u, t);
    valGp(u, mainElement, frontierElement);
    physFluxCu(u, mainElement, frontierElement, flux);
    numFluxUpwind(frontierElement, flux); 
    stiffnessFluxProd(mainElement, flux, SFProd);
    numFluxIntegration(flux, mainElement, frontierElement, fluxVector);
    timeMarching(mainElement, SFProd, fluxVector, step, t, k);

    if(debug) timeChecker(mainElement, frontierElement, flux, u, SFProd, fluxVector, t);
        

}