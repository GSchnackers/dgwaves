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
    u.next.resize(u.node.size(), 0);
    u.numGp.resize(frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    u.bound.resize(mainElement.nodeTags.size(), 0);

    flux.node.resize(mainElement.nodeTags.size() * 3, 0);
    flux.numGp.resize(frontierElement.elementTag.size() * frontierElement.numGp * 3, std::make_pair(0,0));
    flux.direction.resize(flux.numGp.size(), 0);

    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type...");
    setBoundaryConditions(mainElement, u);
    std::cout << "Done." << std::endl;

    SFProd.resize(mainElement.nodeTags.size(), 0);
    fluxVector.resize(mainElement.nodeTags.size(), 0);

    for(t = 0; t < 0.5; t += step)
    {    
        computeBoundaryCondition(mainElement, u, t);
        
        valGp(u, mainElement, frontierElement);
        physFluxCu(u, mainElement, frontierElement, flux);
        numFluxUpwind(frontierElement, flux);
        stiffnessFluxProd(mainElement, flux, SFProd);
        numFluxIntegration(flux, mainElement, frontierElement, fluxVector);
        
        for(i = 0; i < mainElement.elementTag.size(); ++i)
            for(j = 0; j < mainElement.numNodes; ++j)
            {
                int uIndex = i * mainElement.numNodes + j;
                double tmpProd = 0;

                for(k = 0; k < mainElement.numNodes && u.bound[uIndex] > -1; ++k)
                {
                    int matrixIndex = i * mainElement.numNodes * mainElement.numNodes + \
                                 j * mainElement.numNodes + k;
                    int vecIndex = i * mainElement.numNodes + k;

                    tmpProd += mainElement.massMatrixInverse[matrixIndex] * (SFProd[vecIndex] - fluxVector[vecIndex]);
                    
                }

                u.next[uIndex] = u.node[uIndex] + step * tmpProd;

            }
        
        u.node = u.next;

        std::fill(u.next.begin(), u.next.end(), 0);
        std::fill(fluxVector.begin(), fluxVector.end(), 0);
        std::fill(SFProd.begin(), SFProd.end(), 0);

        for(i = 0 ; i < mainElement.elementTag.size(); ++i)
            for(j = 0; j < mainElement.numNodes; ++j)
                mainView.data[i][j] = u.node[i * mainElement.numNodes + j];
            
        gmsh::view::addModelData(mainView.tag, int(t/step), mainView.modelName, mainView.dataType, \
                                 mainElement.elementTag, mainView.data, t, 1);

    }

    
    gmsh::view::write(mainView.tag, "results.msh");

}