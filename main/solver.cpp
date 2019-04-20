// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(Element & mainElement, Element & frontierElement, View & mainView, const double simTime, \
            const double simStep, const int solvType, const int registration, const int debug){

    std::size_t i, j, k; // loop variables.

    double t = 0; // t is the time of the simulation.
    double halfInc = simStep/2;
    double sixthInc = simStep/6;

    Quantity u; // unknowns of the problem.
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

    gmsh::view::addModelData(mainView.tag, int(t/simStep), mainView.modelName, mainView.dataType, \
                                mainElement.elementTag, mainView.data, t, 1);

    gmsh::logger::write("Simulation...");

    // Euler method.
    for(t = 0; t < simTime && !solvType; t += simStep)
    {
        computeCoeff(mainElement, frontierElement, simStep, t, u, flux, k1, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/mainElement.numNodes][i % mainElement.numNodes] = \
            u.node[i] += simStep * k1[i];

        if(!((int(t/simStep) + 1) % registration))
            gmsh::view::addModelData(mainView.tag, int(t/simStep) + 1, mainView.modelName, mainView.dataType, \
                                    mainElement.elementTag, mainView.data, t, 1);

    }

    // Runge-Kutta 4 method.
    for(t = 0; t < simTime && solvType; t += simStep)
    {    
        Quantity uTmp = u;
        computeCoeff(mainElement, frontierElement, simStep, t, uTmp, flux, k1, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k1[i]);
        computeCoeff(mainElement, frontierElement, simStep, t + halfInc, uTmp, flux, k2, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k2[i]);
        computeCoeff(mainElement, frontierElement, simStep, t + halfInc, uTmp, flux, k3, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + simStep * k3[i]);
        computeCoeff(mainElement, frontierElement, simStep, t + simStep, uTmp, flux, k4, debug);

        for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/mainElement.numNodes][i % mainElement.numNodes] = \
            u.node[i] += sixthInc * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]);

        if(!((int(t/simStep) + 1) % registration))
            gmsh::view::addModelData(mainView.tag, int(t/simStep) + 1, mainView.modelName, mainView.dataType, \
                                    mainElement.elementTag, mainView.data, t, 1);

        
    }

    gmsh::logger::write("Done.");
    
    gmsh::view::write(mainView.tag, "results.msh");

}