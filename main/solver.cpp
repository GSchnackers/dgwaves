// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(Element & mainElement, Element & frontierElement, View & mainView, const double simTime, \
            const double simStep, const int solvType, const int registration, const int debug, \
            const double alpha, const std::string & boundFileName, const std::string & propFileName){

    std::size_t i, j, k; // loop variables.

    double t = 0; // t is the time of the simulation.
    double halfInc = simStep/2;
    double sixthInc = simStep/6;

    Quantity u; // unknowns of the problem.
    Quantity flux; // fluxs of the problem.
    Quantity impedance; // impedances of the middle in which each node is.
    Quantity relPermittivity; // Relative permittivities.
    Quantity relPermeability; // Relative permeabilities.
    Quantity conductivity; // Conductivity of the material.

    std::vector<Parameter> bcParam;

    
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
    flux.node.resize(u.node.size() * 3, 0);
    flux.gp.resize(u.gp.size() * 3, std::make_pair(0,0));
    flux.direction.resize(flux.gp.size(), 0);
    flux.num.resize(flux.gp.size(), 0);
    gmsh::logger::write("Done.");

    gmsh::logger::write("Initializing the quantity flux...");
    impedance.gp.resize(mainElement.nodeTags.size());
    gmsh::logger::write("Done.");

    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type...");
    setBoundaryConditions(mainElement, boundFileName, propFileName, relPermittivity, relPermeability, \
                          conductivity, u, bcParam);
    gmsh::logger::write("Done.");

    for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/mainElement.numNodes][i % mainElement.numNodes] = u.node[i];

    gmsh::view::addModelData(mainView.tag, int(t/simStep), mainView.modelName, mainView.dataType, \
                                mainElement.elementTag, mainView.data, t, 1);

    gmsh::logger::write("Simulation...");

    // Euler method.
    for(t = 0; t < simTime && !solvType; t += simStep)
    {
        computeCoeff(mainElement, frontierElement, bcParam, simStep, t, impedance, u, flux, k1, debug, alpha);

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
        computeCoeff(mainElement, frontierElement, bcParam, simStep, t, impedance, uTmp, flux, k1, debug, alpha);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k1[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simStep, t + halfInc, impedance, uTmp, flux, k2,\
                     debug, alpha);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k2[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simStep, t + halfInc, impedance, uTmp, flux, k3, \
                     debug, alpha);

        for (i = 0; i < mainElement.nodeTags.size(); ++i) uTmp.node[i] = u.node[i] * (1 + simStep * k3[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simStep, t + simStep, impedance, uTmp, flux, k4, \
                     debug, alpha);

        for (i = 0; i < mainElement.nodeTags.size(); ++i)
            mainView.data[i/(mainElement.numNodes * 6)][i % (mainElement.numNodes * 6)] = \
            u.node[i] += sixthInc * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]);

        if(!((int(t/simStep) + 1) % registration) || t == simTime - simStep)
            gmsh::view::addModelData(mainView.tag, int(t/simStep) + 1, mainView.modelName, mainView.dataType, \
                                    mainElement.elementTag, mainView.data, t, 6);

        
    }

    gmsh::logger::write("Done.");
    
    gmsh::view::write(mainView.tag, "results.msh");

}