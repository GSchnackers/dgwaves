/* 

    This file initializes the sizes of all the numerical quantities required for solving the PDE with the DG-FEM method.

*/
#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "solver.hpp"
#include "structures.hpp"
#include "parameters.hpp"

void numericalInitializer(const Element & mainElement, Element & frontierElement, \
                          const Simulation & simulation, const PhysicalGroups & physicalGroups,\
                          Quantity & u, Quantity & flux, Properties & matProp){

    std::size_t i;

    gmsh::logger::write("Initializing the BC's parameters...");
    frontierElement.bcParam.resize(9 * frontierElement.elementTag.size());
    gmsh::logger::write("Done.");

    // Initialization of u, which contains the unknowns.
    gmsh::logger::write("Initializing the quantity u...");
    u.node.resize(simulation.uNum * mainElement.nodeTags.size(), 0);
    u.gp.resize(simulation.uNum * frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    gmsh::logger::write("Done.");


    gmsh::logger::write("Initializing the quantity flux...");
    flux.node.resize(3 * simulation.uNum * mainElement.nodeTags.size(), 0);
    flux.gp.resize(3 * simulation.uNum * frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    flux.direction.resize(flux.gp.size(), 0);
    flux.num.resize(flux.gp.size()/3, 0);
    gmsh::logger::write("Done.");

    gmsh::logger::write("Initializing the quantity impedance and inductances...");
    matProp.speedGpInv.resize(frontierElement.numGp * frontierElement.elementTag.size());
    matProp.speedGpSumInvInv.resize(matProp.speedGpInv.size());
    gmsh::logger::write("Done.");

    // Loading the physical properties of the material.
    gmsh::logger::write("Initializing the quantity conductivity...");
    matProp.conductivity.node.resize(mainElement.nodeTags.size(), 0);
    matProp.conductivity.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, \
                                   std::make_pair(0,0));
    gmsh::logger::write("Done.");

    // Loading the physical properties of the material.
    gmsh::logger::write("Initializing the quantity relative permittivity...");
    matProp.relPermittivity.node.resize(mainElement.nodeTags.size(), 0);
    matProp.relPermittivity.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, \
                                      std::make_pair(0,0));
    gmsh::logger::write("Done.");

    // Loading the physical properties of the material.
    gmsh::logger::write("Initializing the quantity relative permeability...");
    matProp.relPermeability.node.resize(mainElement.nodeTags.size(), 0);
    matProp.relPermeability.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, \
                                      std::make_pair(0,0));
    gmsh::logger::write("Done.");
    
    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type and the material properties...");
    setProperties(mainElement, frontierElement, simulation, physicalGroups, matProp);
    setBoundaryCondition(frontierElement, simulation, physicalGroups, u);
    gmsh::logger::write("Done.");

    
    // Compute the value of the material properties at the Gauss points.
    gmsh::logger::write("Computing the relative permeability at the Gauss points.");
    valGp(matProp.relPermeability, mainElement, frontierElement, 1, matProp);
    gmsh::logger::write("Done.");
    gmsh::logger::write("Computing the relative permittivity at the Gauss points.");
    valGp(matProp.relPermittivity, mainElement, frontierElement, 1, matProp);
    gmsh::logger::write("Done.");

    // Computes the adimensionnal coefficients.
    gmsh::logger::write("Setting the adimensionnal numbers, the impedance and conductances at the Gauss points\
                         and at the nodes...");
    for(i = 0; i < matProp.speedGpInv.size(); ++i)
    {   
        matProp.speedGpInv[i].first = sqrt(matProp.relPermittivity.gp[i].first * matProp.relPermeability.gp[i].first);
        matProp.speedGpInv[i].second = sqrt(matProp.relPermittivity.gp[i].second * matProp.relPermeability.gp[i].second);
        matProp.speedGpSumInvInv[i] = 1/(matProp.speedGpInv[i].first + matProp.speedGpInv[i].second);       
    }

    gmsh::logger::write("Done");

}