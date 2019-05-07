#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void numericalInitializer(const Element & mainElement, Element & frontierElement, \
                          const Simulation & simulation, const PhysicalGroups & physicalGroups,\
                          Quantity & u, Quantity & flux, Properties & matProp, \
                          std::vector<Parameter> & bcParam){

    std::size_t i;

    double vacuumImp = std::sqrt(VACUUM_PERMEA/VACUUM_PERMIT);

    gmsh::logger::write("Initializing the bc parameters...");
    bcParam.resize(simulation.uNum * mainElement.nodeTags.size());
    for(i = 0; i < bcParam.size(); ++i)
    {
        bcParam[i].param1 = 0;
        bcParam[i].param2 = 0;
        bcParam[i].param3 = 0;
    }
    gmsh::logger::write("Done.");

    // Initialization of u, which contains the unknowns.
    gmsh::logger::write("Initializing the quantity u...");
    u.node.resize(simulation.uNum * mainElement.nodeTags.size(), 0);
    u.gp.resize(simulation.uNum * frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    u.bound.resize(simulation.uNum * mainElement.nodeTags.size(), 0);
    u.boundSign.resize(simulation.uNum * mainElement.nodeTags.size(), 0);
    gmsh::logger::write("Done.");


    gmsh::logger::write("Initializing the quantity flux...");
    flux.node.resize(3 * simulation.uNum * mainElement.nodeTags.size(), 0);
    flux.gp.resize(3 * simulation.uNum * frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    flux.direction.resize(flux.gp.size(), 0);
    flux.num.resize(flux.gp.size(), 0);
    gmsh::logger::write("Done.");

    gmsh::logger::write("Initializing the quantity impedance and inductances...");
    matProp.impedance.node.resize(mainElement.nodeTags.size(), 0);
    matProp.impedance.bound.resize(mainElement.nodeTags.size(), 0);
    matProp.impedance.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    matProp.conductance.node.resize(mainElement.nodeTags.size(), 0);
    matProp.conductance.bound.resize(mainElement.nodeTags.size(), 0);
    matProp.conductance.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    gmsh::logger::write("Done.");

    // Loading the physical properties of the material.
    gmsh::logger::write("Initializing the quantity conductivity...");
    matProp.conductivity.node.resize(mainElement.nodeTags.size(), 0);
    matProp.conductivity.bound.resize(mainElement.nodeTags.size(), 0);
    matProp.conductivity.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, \
                                   std::make_pair(0,0));
    gmsh::logger::write("Done.");

    // Loading the physical properties of the material.
    gmsh::logger::write("Initializing the quantity relative permittivity...");
    matProp.relPermittivity.node.resize(mainElement.nodeTags.size(), 0);
    matProp.relPermittivity.bound.resize(mainElement.nodeTags.size(), 0);
    matProp.relPermittivity.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, \
                                      std::make_pair(0,0));
    gmsh::logger::write("Done.");

    // Loading the physical properties of the material.
    gmsh::logger::write("Initializing the quantity relative permeability...");
    matProp.relPermeability.node.resize(mainElement.nodeTags.size(), 0);
    matProp.relPermeability.bound.resize(mainElement.nodeTags.size(), 0);
    matProp.relPermeability.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, \
                                      std::make_pair(0,0));
    gmsh::logger::write("Done.");

    gmsh::logger::write("Initializing the adimensionnal quantity...");
    matProp.eta.node.resize(mainElement.nodeTags.size(), 0);
    matProp.eta.bound.resize(mainElement.nodeTags.size(), 0);
    matProp.eta.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp, std::make_pair(0,0));
    gmsh::logger::write("Done.");

    // Setting of the properties of the elements
    setProperties(mainElement, frontierElement, simulation, physicalGroups, matProp);

    
    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type and the material properties...");
    setBoundaryCondition(frontierElement, mainElement, simulation, physicalGroups, u, bcParam);
    gmsh::logger::write("Done.");

    
    // Compute the value of the material properties at the Gauss points.
    gmsh::logger::write("Computing the relative permeability at the Gauss points.");
    valGp(matProp.relPermeability, mainElement, frontierElement, 1);
    gmsh::logger::write("Done.");
    gmsh::logger::write("Computing the relative permittivity at the Gauss points.");
    valGp(matProp.relPermittivity, mainElement, frontierElement, 1);
    gmsh::logger::write("Done.");

    // Computes the adimensionnal coefficients.
    gmsh::logger::write("Setting the adimensionnal numbers, the impedance and conductances at the Gauss points\
                         and at the nodes...");
    for(i = 0; i < mainElement.nodeTags.size(); ++i)
    {
        
        matProp.impedance.bound[i] = matProp.impedance.node[i] = \
                vacuumImp * std::sqrt(matProp.relPermittivity.node[i]/matProp.relPermeability.node[i]);

        matProp.conductance.bound[i] = matProp.conductance.node[i] = 1/matProp.impedance.node[i];

        matProp.eta.bound[i] = matProp.eta.node[i] = simulation.L * matProp.conductivity.node[i] * \
                              matProp.impedance.node[i]/matProp.relPermittivity.node[i];
                              
    }

    // Value at the Gauss points of the various quantities.
    valGp(matProp.impedance, mainElement, frontierElement, 1);
    valGp(matProp.conductance, mainElement, frontierElement, 1);
    valGp(matProp.eta, mainElement, frontierElement, 1);

    gmsh::logger::write("Done");

}