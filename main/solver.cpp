// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void solver(const Element & mainElement, const Element & frontierElement, const PhysicalGroups & physicalGroups,\
            View & mainView, Simulation & simulation){

    std::size_t i, j, k; // loop variables.

    double t = 0; // t is the time of the simulation.
    double halfInc = simulation.simStep/2;
    double sixthInc = simulation.simStep/6;

    Quantity u; // unknowns of the problem.
    Quantity flux; // fluxs of the problem.
    Simulation simParam; // Parameters of the simulation.
    Properties matProp; // Properties of the material all over the domain.
    std::vector<Parameter> bcParam; // Parameters associated to the boundary conditions.

    std::vector<double> k1(6 * mainElement.nodeTags.size(), 0);
    std::vector<double> k2(6 * mainElement.nodeTags.size(), 0);
    std::vector<double> k3(6 * mainElement.nodeTags.size(), 0);
    std::vector<double> k4(6 * mainElement.nodeTags.size(), 0);

    // Initializes the size of the useful quantities.
    numericalInitializer(mainElement, frontierElement, simulation, physicalGroups, u, flux, matProp, bcParam);
    
    for (i = 0; i < u.node.size(); ++i)
        mainView.data[i/(6 * mainElement.numNodes)][i % (6 * mainElement.numNodes)] = u.node[i];

    std::cout << "Hello" << std::endl;

    gmsh::view::addModelData(mainView.tag, int(t/simulation.simStep), mainView.modelName, mainView.dataType, \
                                mainElement.elementTag, mainView.data, t, 6);


    gmsh::logger::write("Simulation...");

    // Euler method.
    for(t = 0; t < simulation.simTime && !simulation.solver; t += simulation.simStep)
    {
        
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t, u, flux, k1);

        for (i = 0; i < u.node.size(); ++i)
            mainView.data[i/(mainElement.numNodes * 6)][i % (mainElement.numNodes * 6)] = \
                                                  u.node[i] += simulation.simStep * k1[i];

        if(!((int(t/simulation.simStep) + 1) % simulation.registration))
            gmsh::view::addModelData(mainView.tag, int(t/simulation.simStep) + 1, mainView.modelName, \
                                    mainView.dataType, mainElement.elementTag, mainView.data, t, 6);

        if(!(int(t/simulation.simTime) % 20))
            std::cout << "\33\rProgression: " << t/simulation.simTime * 100 << "%";

    }

    // Runge-Kutta 4 method.
    for(t = 0; t < simulation.simTime && simulation.solver; t += simulation.simStep)
    {    

        int stepNum = int(t/simulation.simTime);
        Quantity uTmp = u;
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t, uTmp, flux, k1);

        for (i = 0; i < u.node.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k1[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t + halfInc, uTmp, flux, k2);

        for (i = 0; i < u.node.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k2[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t + halfInc, uTmp, flux, k3);

        for (i = 0; i < u.node.size(); ++i) uTmp.node[i] = u.node[i] * (1 + simulation.simStep * k3[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t + simulation.simStep, uTmp,\
                     flux, k4);

        for (i = 0; i < u.node.size(); ++i)
            mainView.data[i/(mainElement.numNodes * 6)][i % (mainElement.numNodes * 6)] = \
            u.node[i] += sixthInc * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]);

        if(!((stepNum + 1) % simulation.registration) || t == simulation.simTime - simulation.simStep)
            gmsh::view::addModelData(mainView.tag, int(t/simulation.simStep) + 1, mainView.modelName,\
                                    mainView.dataType, mainElement.elementTag, mainView.data, t, 6);

        if(!(stepNum % 20))
            std::cout << "\33\rProgression: " << t/simulation.simTime * 100 << "%";
        
    }

    std::cout << std::endl;

    gmsh::logger::write("Done.");
    
    gmsh::view::write(mainView.tag, "results.msh");

}