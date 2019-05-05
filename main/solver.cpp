// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Functions that allows to display the results.
void resultsDisp(const Element & mainElement, const Quantity & u, const Simulation & simulation, double t, \
                 View & EView, View & HView){

    std::size_t i, j, k;

    for (i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            for(k = 0; k < 3; ++k)
            {
                int uEIndex = i * mainElement.numNodes * 6  + j * 6 + k;
                int uHINdex = uEIndex + 3;

                EView.data[i][j * 3 + k] = u.node[uEIndex];
                HView.data[i][j * 3 + k] = u.node[uHINdex];
            }
        
    gmsh::view::addModelData(EView.tag, int(t/simulation.simStep), EView.modelName, EView.dataType, \
                        mainElement.elementTag, EView.data, t, 3);
    gmsh::view::addModelData(HView.tag, int(t/simulation.simStep), HView.modelName, HView.dataType, \
                        mainElement.elementTag, HView.data, t, 3);

}

void solver(const Element & mainElement, const Element & frontierElement, const PhysicalGroups & physicalGroups,\
            View & EView, View & HView, Simulation & simulation){

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
    
    resultsDisp(mainElement, u, simulation, t, EView, HView);

    gmsh::logger::write("Simulation...");

    // Euler method.
    for(t = 0; t < simulation.simTime && !simulation.solver; t += simulation.simStep)
    {
        
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t, u, flux, k1);

        for (i = 0; i < u.node.size(); ++i) u.node[i] += simulation.simStep * k1[i];

        if(!((int(t/simulation.simStep) + 1) % simulation.registration))
            resultsDisp(mainElement, u, simulation, t, EView, HView);

        if(!(int(t/simulation.simTime) % 20) || t/simulation.simTime == 1)
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

        for(i = 0; i < u.node.size(); ++i) u.node[i] += sixthInc * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);

        if(!((int(t/simulation.simStep) + 1) % simulation.registration))
            resultsDisp(mainElement, u, simulation, t, EView, HView);

        if(!(stepNum % 20))
            std::cout << "\33\rProgression: " << t/simulation.simTime * 100 << "%";
        
    }

    std::cout << std::endl;

    gmsh::logger::write("Done.");
    
    gmsh::view::write(EView.tag, "electricalField.msh");
    gmsh::view::write(EView.tag, "magneticField.msh");

}