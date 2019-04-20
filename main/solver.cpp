// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(Element & mainElement, Element & frontierElement, View & mainView){

    std::size_t i, j, k; // loop variables.

    double t, increment = 0.001; // t is the time, increment is the time increment of the simulation.
    double halfInc = increment/2;
    double sixthInc = increment/6;

    Quantity u; // unknowns of the problem.
    Quantity uTmp;
    Quantity flux; // fluxs of the problem.
    
    std::vector<double> k1(mainElement.nodeTags.size(), 0);
    std::vector<double> k2(mainElement.nodeTags.size(), 0);
    std::vector<double> k3(mainElement.nodeTags.size(), 0);
    std::vector<double> k4(mainElement.nodeTags.size(), 0);

    // Initialization of the nodal values.
    gmsh::logger::write("Initializing the quantity u...");
    u.node.resize(mainElement.nodeTags.size(), 0);
    u.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    u.bound.resize(mainElement.nodeTags.size(), 0);
    u.boundSign.resize(mainElement.nodeTags.size(), 0);
    gmsh::logger::write("Done.");

    gmsh::logger::write("Initializing the quantity flux...");
    flux.node.resize(mainElement.nodeTags.size() * 3, 0);
    flux.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp * 3, std::make_pair(0,0));
    flux.direction.resize(flux.gp.size(), 0);
    flux.num.resize(flux.gp.size(), 0);
    flux.bound.resize(mainElement.nodeTags.size() * 3, 0);
    gmsh::logger::write("Done.");

    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type...");
    setBoundaryConditions(mainElement, u);
    gmsh::logger::write("Done.");

    for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/mainElement.numNodes][i % mainElement.numNodes] = u.node[i];

    gmsh::view::addModelData(mainView.tag, int(t/increment), mainView.modelName, mainView.dataType, \
                                mainElement.elementTag, mainView.data, t, 1);

    gmsh::logger::write("Simulation...");
    for(t = 0; t < .5; t += increment)
    {    
        uTmp = u;
        computeCoeff(mainElement, frontierElement, increment, t, uTmp, flux, k1);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] *(1 + halfInc * k1[i]);
        computeCoeff(mainElement, frontierElement, increment, t + halfInc, uTmp, flux, k2);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] *(1 + halfInc * k2[i]);
        computeCoeff(mainElement, frontierElement, increment, t + halfInc, uTmp, flux, k3);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] *(1 + increment * k3[i]);
        computeCoeff(mainElement, frontierElement, increment, t + increment, uTmp, flux, k4);

        for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/mainElement.numNodes][i % mainElement.numNodes] = \
            u.node[i] += sixthInc * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]);

        gmsh::view::addModelData(mainView.tag, int(t/increment) + 1, mainView.modelName, mainView.dataType, \
                                mainElement.elementTag, mainView.data, t, 1);
        
    }
    
    gmsh::logger::write("Done.");
    
    gmsh::view::write(mainView.tag, "results.msh");

}