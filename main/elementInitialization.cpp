/*

    This functions initializes the caracteristics of an element of type elementType and of dimension
    dim.

*/

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void elementInitialization(struct Element & element2D){

    gmsh::model::mesh::getElementTypes(element2D.elementTypes, 2); // Gets the type of all elements of dimension 2 of all tags
    if(element2D.elementTypes.size() > 1){
        gmsh::logger::write("No hybrid implementation is authorized.", "error");
        exit(-1);
    }

    gmsh::model::mesh::getElementProperties(element2D.elementTypes[0], element2D.name, \
                                            element2D.dim, element2D.order, element2D.numNodes,\
                                            element2D.parametricCoordinates); // Load the element properties in the memory.

}

// REM: can be easily modified to comprehend hybrid elements.