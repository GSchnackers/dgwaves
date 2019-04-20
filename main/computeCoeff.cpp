#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void computeCoeff(const Element & mainElement, const Element & frontierElement, const double step, \
                  const double t, Quantity & u, Quantity & flux, std::vector<double> & k){

    std::vector<double> SFProd(mainElement.nodeTags.size(), 0);
    std::vector<double> fluxVector(mainElement.nodeTags.size(), 0);

    computeBoundaryCondition(mainElement, u, t);
    valGp(u, mainElement, frontierElement);
    physFluxCu(u, mainElement, frontierElement, flux);
    numFluxUpwind(frontierElement, flux); 
    stiffnessFluxProd(mainElement, flux, SFProd);
    numFluxIntegration(flux, mainElement, frontierElement, fluxVector);
    timeMarching(mainElement, SFProd, fluxVector, step, t, k);

}