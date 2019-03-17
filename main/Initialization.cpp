#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Function that initializes the characteristics of the elements represented by element.

void Initialization(Element & element, const int meshDim, const int gaussType){

    std::size_t i, j; // Index variables

    gmsh::vectorpair entities; // Vector pair used to contain the dimension and tags of entities.

    std::vector<double> bin1, bin2; // Variable used to collect useless results.

    // Integration Type definition.
    std::string integrationType = "Gauss" + std::to_string(gaussType);

    // Loading of the entities.
    gmsh::model::getEntities(entities, meshDim); // Load the entities of dimension meshDim tag in the memory.

    // Initialization of the 2D elements.

    gmsh::model::mesh::getElementTypes(element.elementType, meshDim);
    if(element.elementType.size() > 1){

        gmsh::logger::write("No hybrid implementation is authorized.\
                                The program will now exit.", "error");
        exit(-1);

    }

    // Gets the elements tag of the entities, the node tags of the geometry.
    gmsh::model::mesh::getElementsByType(element.elementType[0], element.elementTag, element.nodeTags);

    // Gets the frontier nodes of the main elements. Depends on the mesh dimension. 
    if(meshDim == 3)    
        gmsh::model::mesh::getElementFaceNodes(element.elementType[0], 3, element.frontierNode); 
    else if(meshDim == 2)                             
        gmsh::model::mesh::getElementEdgeNodes(element.elementType[0], element.frontierNode);
    else if(meshDim == 1)
        gmsh::model::mesh::getElementEdgeNodes(element.elementType[0], element.frontierNode, -1, true);

    // Gets the number of nodes per edge of an elements.

    for(j = 0; j < element.frontierNode.size(); ++j)
        if(element.frontierNode[j] == element.frontierNode[0]){

            element.numberEdgeNode = j + 1;
            break;

        }

    // Gets the shape functions of the 2D elements at the Gauss points;
    gmsh::model::mesh::getBasisFunctions(element.elementType[0], integrationType, "IsoParametric",\
                                element.gaussPointsParam, element.numComponentShapeFunctions,\
                                element.shapeFunctionsParam);

    // Gets the gradient of the shape functions of the 2D elements at the Gauss points;
    gmsh::model::mesh::getBasisFunctions(element.elementType[0], integrationType, "GradIsoParametric",\
                                            bin1, element.numComponentShapeFunctionsGrad,\
                                            element.shapeFunctionsGradParam);

    // Gets the jacobian at each Gauss point of the 2D elements.
    gmsh::model::mesh::getJacobians(element.elementType[0], integrationType, bin1, element.jacobians,\
                                    bin2);

    // Gets the type of all elements of dimension 2 of all tags.

}