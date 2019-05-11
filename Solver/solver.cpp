// This file contains the solver of the problem.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "solver.hpp"
#include "structures.hpp"

void sinusoidalDisp(const Element & mainElement, const Quantity & u, const Simulation & simulation, double t, \
                 View & view){

    std::size_t i, j;

    #pragma omp parallel for shared(u, view, i) private(j)
    for(i = 0; i < mainElement.elementTag.size(); ++i)
        for(j = 0; j < mainElement.numNodes; ++j)
            view.data[i][j] = u.node[i * mainElement.numNodes + j];
    
    gmsh::view::addModelData(view.tag, int(t/simulation.simStep), view.modelName, view.dataType, \
                        mainElement.elementTag, view.data, t, 1);

}

// Functions that allows to display the results.
void ELMDisp(const Element & mainElement, const Quantity & u, const Simulation & simulation, double t, \
                 View & EView, View & HView){

    std::size_t i, j, k;

    for(i = 0; i < mainElement.elementTag.size(); ++i)
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

void solver(const Element & mainElement, Element & frontierElement, const PhysicalGroups & physicalGroups,\
            View & view1, View & view2, Simulation & simulation){

    std::size_t i, j, k; // loop variables.

    double t = 0; // t is the time of the simulation.
    double halfInc = simulation.simStep/2;
    double sixthInc = simulation.simStep/6;

    Quantity u; // unknowns of the problem.
    Quantity flux; // fluxs of the problem.
    Simulation simParam; // Parameters of the simulation.
    Properties matProp; // Properties of the material all over the domain.
    std::vector<Parameter> bcParam; // Parameters associated to the boundary conditions.
    
    std::vector<double> error(simulation.uNum * int(simulation.simTime / simulation.simStep), 0);
    std::vector<double> errorNodes(simulation.uNum * mainElement.nodeTags.size(), 0);
    std::vector<double> coordinates(3 * mainElement.nodeTags.size(), 0);
    std::vector<double> coordNodes(3 * mainElement.numNodes, 0);
    std::vector<double> binParam(3 * mainElement.numNodes);

    std::vector<double> k1(simulation.uNum * mainElement.nodeTags.size(), 0);
    std::vector<double> k2(simulation.uNum * mainElement.nodeTags.size(), 0);
    std::vector<double> k3(simulation.uNum * mainElement.nodeTags.size(), 0);
    std::vector<double> k4(simulation.uNum * mainElement.nodeTags.size(), 0);

    for(i = 0; i < mainElement.nodeTags.size(); i++){
        gmsh::model::mesh::getNode(mainElement.nodeTags[i], coordNodes, binParam);

        for(j = 0; j < 3; j++){
            coordinates[3*i + j] = coordNodes[j];
        }
    }

    // Initializes the size of the useful quantities.
    numericalInitializer(mainElement, frontierElement, simulation, physicalGroups, u, flux, matProp, bcParam);
    
    if(simulation.uNum == 6)
        ELMDisp(mainElement, u, simulation, t, view1, view2);

    else if(simulation.uNum == 1)
        sinusoidalDisp(mainElement, u, simulation, t, view1);

    gmsh::logger::write("Simulation...");

    // Euler method.
    for(t = 0; t < simulation.simTime && !simulation.solver; t += simulation.simStep)
    {
        
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t, u, flux, k1);

        #pragma omp parallel for shared(u, i, k1)
        for (i = 0; i < u.node.size(); ++i) u.node[i] += simulation.simStep * k1[i];

        if(simulation.error){
            compare(error[int(t/simulation.simStep)], errorNodes, u, coordinates, mainElement, simulation, bcParam, t);
        }

        if(!((int(t/simulation.simStep) + 1) % simulation.registration))
        {
            if(simulation.uNum == 6)
                ELMDisp(mainElement, u, simulation, t, view1, view2);

            else if(simulation.uNum == 1)
                sinusoidalDisp(mainElement, u, simulation, t, view1);
        }

        if(!(int(t/simulation.simTime) % 20) || t/simulation.simTime == 1)
            std::cout << "\33\rProgression: " << int(t/simulation.simTime * 100) << "%";

    }

    // Runge-Kutta 4 method.
    for(t = 0; t < simulation.simTime && simulation.solver; t += simulation.simStep)
    {    

        int stepNum = int(t/simulation.simTime);
        Quantity uTmp = u;
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t, uTmp, flux, k1);

        #pragma omp parallel for shared(u, i, k1, uTmp)
        for (i = 0; i < u.node.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k1[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t + halfInc, uTmp, flux, k2);

        #pragma omp parallel for shared(u, i, k2, uTmp)
        for (i = 0; i < u.node.size(); ++i) uTmp.node[i] = u.node[i] * (1 + halfInc * k2[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t + halfInc, uTmp, flux, k3);

        #pragma omp parallel for shared(u, i, k3, uTmp)
        for (i = 0; i < u.node.size(); ++i) uTmp.node[i] = u.node[i] * (1 + simulation.simStep * k3[i]);
        computeCoeff(mainElement, frontierElement, bcParam, simulation, matProp, t + simulation.simStep, uTmp,\
                     flux, k4);

        #pragma omp parallel for shared(u, i, k1, k2, k3, k4, uTmp)
        for(i = 0; i < u.node.size(); ++i) u.node[i] += sixthInc * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);

        if(simulation.error)
            compare(error[int(t/simulation.simStep)], errorNodes, u, coordinates, mainElement, simulation, bcParam, t);
        
        if(!((int(t/simulation.simStep) + 1) % simulation.registration))
        {
            if(simulation.uNum == 6)
                ELMDisp(mainElement, u, simulation, t, view1, view2);

            else if(simulation.uNum == 1)
                sinusoidalDisp(mainElement, u, simulation, t, view1);
        }

        if(!(stepNum % 20))
            std::cout << "\33\rProgression: " << int(t/simulation.simTime * 100) << "%";
        
    }

    std::cout << std::endl;

    gmsh::logger::write("Done.");
    
    if(simulation.error)
        writeError(error, simulation);

    gmsh::view::write(view1.tag, "electricField.msh");

    if(simulation.uNum == 6)
        gmsh::view::write(view2.tag, "magneticField.msh");

}