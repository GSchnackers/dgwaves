// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(Element & mainElement, Element & frontierElement, View & mainView, const double simTime, \
            const double incrementation, const int solvType, const int registration, const int debug){

    std::size_t i, j, k; // loop variables.

    double t, increment = 0.001; // t is the time, increment is the time increment of the simulation.
    double halfInc = increment/2;
    double sixthInc = increment/6;

    Quantity u; // unknowns of the problem.
    Quantity flux; // fluxs of the problem.

    
    std::vector<double> k1(6 * mainElement.nodeTags.size(), 0);
    std::vector<double> k2(6 * mainElement.nodeTags.size(), 0);
    std::vector<double> k3(6 * mainElement.nodeTags.size(), 0);
    std::vector<double> k4(6 * mainElement.nodeTags.size(), 0);

    // Initialization of the nodal values.
    gmsh::logger::write("Initializing the quantity u...");
    u.node.resize(6 * mainElement.nodeTags.size(), 0);
    u.gp.resize(6 * frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    u.bound.resize(6 * mainElement.nodeTags.size(), 0);
    u.boundSign.resize(6 * mainElement.nodeTags.size(), 0);
    gmsh::logger::write("Done.");

    gmsh::logger::write("Initializing the quantity flux...");
    flux.node.resize(mainElement.nodeTags.size() * 18, 0);
    flux.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp * 18, std::make_pair(0,0));
    flux.direction.resize(flux.gp.size(), 0);
    flux.num.resize(flux.gp.size(), 0);
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

    // Euler method.
    for(t = 0; t < simTime && !solvType; t += increment)
    {
        computeCoeff(mainElement, frontierElement, increment, t, u, flux, k1, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/mainElement.numNodes][i % mainElement.numNodes] = \
            u.node[i] += increment * k1[i];

        if(!((int(t/increment) + 1) % registration))
            gmsh::view::addModelData(mainView.tag, int(t/increment) + 1, mainView.modelName, mainView.dataType, \
                                    mainElement.elementTag, mainView.data, t, 1);

    }

    // Runge-Kutta 4 method.
    for(t = 0; t < simTime && solvType; t += increment)
    {    
        Quantity uTmp = u;
        computeCoeff(mainElement, frontierElement, increment, t, uTmp, flux, k1, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k1[i]);
        computeCoeff(mainElement, frontierElement, increment, t + halfInc, uTmp, flux, k2, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k2[i]);
        computeCoeff(mainElement, frontierElement, increment, t + halfInc, uTmp, flux, k3, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + increment * k3[i]);
        computeCoeff(mainElement, frontierElement, increment, t + increment, uTmp, flux, k4, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/mainElement.numNodes][i % mainElement.numNodes] = \
            u.node[i] += sixthInc * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]);

        if(!((int(t/increment) + 1) % registration) || t == simTime - 1)
            gmsh::view::addModelData(mainView.tag, int(t/increment) + 1, mainView.modelName, mainView.dataType, \
                                    mainElement.elementTag, mainView.data, t, 1);

        
    }

    gmsh::logger::write("Done.");
    
    gmsh::view::write(mainView.tag, "results.msh");

}