// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(Element & mainElement, Element & frontierElement){

    double t; // t is the time.
    std::size_t i, j; // loop variables.

    Quantity u; // unknowns of the problem.
    Quantity flux; // fluxs of the problem.

    std::vector<double> SFProd(mainElement.nodeTags.size(), 0);

    // Initialization of the nodal values.
    u.node.resize(mainElement.nodeTags.size());

    for(t = 0; t = 1; t += 0.01){

        computeBoundaryCondition(frontierElement, u, t);
        physFluxCu(u, mainElement, frontierElement, flux);
        numFluxUpwind(frontierElement, flux);
        stiffnessFluxProd(mainElement, flux, SFProd);

    }

}