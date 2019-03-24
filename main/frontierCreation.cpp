#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void frontierCreation(const Element mainElement, Element & frontierElement, const int meshDim,\
                      const std::vector<int> sortedNodes){

    // Allocation for element type.
    frontierElement.elementType.resize(1);
    
    // Addition of a new entity to the model of dimension meshDim - 1, which is the dimension of the frontiers of the domain.
    frontierElement.entityTag = gmsh::model::addDiscreteEntity(meshDim - 1);
    std::string frontierName;

    // Selection of the name of the frontier.
    if(meshDim == 3) frontierName = "triangle";
    else if(meshDim == 2) frontierName = "line";

    // Gets the type of the frontier elements.
    frontierElement.elementType[0] = gmsh::model::mesh::getElementType(frontierName, mainElement.order);

    // Set the elements in the model associated to the new discrete entity.
    gmsh::model::mesh::setElementsByType(meshDim - 1, frontierElement.entityTag, frontierElement.elementType[0], {},\
                                         sortedNodes);

}