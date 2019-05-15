#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "solver.hpp"
#include "structures.hpp"
#include "parameters.hpp"

#define VACUUM_PERMIT 8.54e-12
#define VACUUM_PERMEA 1.257e-6
#define VACUUM_IMPEDANCE 377
#define VACUUM_CONDUCTANCE 2.652e-3

void numericalInitializer(const Element & mainElement, Element & frontierElement, \
                          const Simulation & simulation, const PhysicalGroups & physicalGroups,\
                          Quantity & u, Quantity & flux, Properties & matProp){

    std::size_t i;

    double vacuumImp = std::sqrt(VACUUM_PERMEA/VACUUM_PERMIT);

    gmsh::logger::write("Initializing the bc parameters...");
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
    matProp.speedGp.resize(frontierElement.numGp * frontierElement.elementTag.size());
    matProp.speedGpSumInv.resize(matProp.speedGp.size());
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

    gmsh::logger::write("Initializing the adimensionnal quantity...");
    matProp.eta.node.resize(mainElement.nodeTags.size(), 0);
    matProp.eta.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    gmsh::logger::write("Done.");

    // Setting of the properties of the elements
    setProperties(mainElement, frontierElement, simulation, physicalGroups, matProp);

    
    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type and the material properties...");
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
    for(i = 0; i < matProp.speedGp.size(); ++i)
    {   
        matProp.speedGp[i].first = sqrt(matProp.relPermittivity.gp[i].first * matProp.relPermeability.gp[i].first);
        matProp.speedGp[i].second = sqrt(matProp.relPermittivity.gp[i].second * matProp.relPermeability.gp[i].second);
        matProp.speedGpSumInv[i] = 1/(matProp.speedGp[i].first + matProp.speedGp[i].second);       
    }

    gmsh::logger::write("Done");

}