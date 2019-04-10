// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(Element & mainElement, Element & frontierElement, View & mainView){

    double t, step = 0.01; // t is the time.
    std::size_t i, j, k; // loop variables.

    Quantity u; // unknowns of the problem.
    Quantity flux; // fluxs of the problem.

    std::vector<double> SFProd;
    std::vector<double> fVector;

    std::vector<std::vector<double>> data(mainElement.nodeTags.size(), std::vector<double>(1));

    // Initialization of the nodal values.
    u.node.resize(mainElement.nodeTags.size());
    u.next.resize(u.node.size());
    u.bound.resize(frontierElement.nodeTags.size());
    u.numGp.resize(frontierElement.elementTag.size() * frontierElement.numGp);

    flux.node.resize(mainElement.nodeTags.size() * 3);
    flux.numGp.resize(frontierElement.elementTag.size() * frontierElement.numGp * 3);

    SFProd.resize(mainElement.nodeTags.size());
    fVector.resize(mainElement.nodeTags.size());

    for(t = 0; t < 25 * step; t += step){

        computeBoundaryCondition(mainElement, frontierElement, u, t);
        valGp(u, mainElement, frontierElement);
        physFluxCu(u, mainElement, frontierElement, flux);
        numFluxUpwind(frontierElement, flux);
        stiffnessFluxProd(mainElement, flux, SFProd);
        numFluxIntegration(flux, mainElement, frontierElement, fVector);

        for(i = 0; i < mainElement.elementTag.size(); ++i)
            for(j = 0; j < mainElement.numNodes; ++j)
            {
                int uIndex = i * mainElement.numNodes + j;
                u.next[uIndex] = 0;

                for(k = 0; k < mainElement.numNodes; ++k)
                {
                    int mIndex = i * mainElement.numNodes * mainElement.numNodes + j * mainElement.numNodes + k;
                    int vecIndex = i * mainElement.numNodes + k;

                    u.next[uIndex] += u.node[uIndex] + step * mainElement.massMatrixInverse[mIndex] * \
                                     (SFProd[vecIndex] + fVector[vecIndex]);


                }


            }
        
        u.node = u.next;

        for(i = 0 ; i < u.node.size(); ++i)
            data[i][0] = u.node[i];

        std::cout << u.node.size() << std::endl;

        gmsh::view::addModelData(mainView.tag, int(t/0.01), mainView.modelName, mainView.dataType, \
                                 mainElement.nodeTags, data, t, 1);
        gmsh::view::write(mainView.tag, "results.msh");

    }

}