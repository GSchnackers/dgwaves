// Header that deals with the parameters functions.

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "structures.hpp"

/*
   Function that reads the parameters of the simulation. The parameters are read in a .wave file.
*/
void readParam(std::string fileName, Simulation & simulation);

/*
   Function that set the properties of the different materials on the domain. It is udes with a .prop file.
*/
void setProperties(const Element & mainElement, const Element & frontierElement, const Simulation & simulation, \
                   const PhysicalGroups & physicalGroups, Properties & matProp);

#endif