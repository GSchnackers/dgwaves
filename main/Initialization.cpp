#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Function that initializes the characteristics of the elements represented by element.

void Initialization(Element & element, const int meshDim, const int gaussType){

    std::size_t i, j; // Index variables

    gmsh::vectorpair entities; // Vector pair used to contain the dimension and tags of entities.

    std::vector<double> bin1, bin2, bin3; // Variable used to collect useless results.

    // Integration Type definition.
    std::string integrationType = "Gauss" + std::to_string(gaussType);

    // Loading of the entities.
    gmsh::model::getEntities(entities, meshDim); // Load the entities of dimension meshDim tag in the memory.

    // Initialization of the elements.
    gmsh::model::mesh::getElementTypes(element.elementType, meshDim);
    if(element.elementType.size() > 1){

        gmsh::logger::write("No hybrid implementation is authorized.\
                                The program will now exit.", "error");
        exit(-1);

    }

    // Gets the properties of the element.
    gmsh::model::mesh::getElementProperties(element.elementType[0], element.name, element.dim,\
                                            element.order, element.numNodes,\
                                            element.parametricCoordinates);

    // Get the number of sides of the element from its name. Its name is the family name + a number.
    if(element.name.find("Triangle") != std::string::npos) element.numSide = 3;
    else if(element.name.find("Line") != std::string::npos) element.numSide = 1;
    else if(element.name.find("Quadrangle") != std::string::npos) element.numSide = 4;
    else if(element.name.find("Tetrahedron") != std::string::npos) element.numSide = 6;
    else if(element.name.find("Pyramid") != std::string::npos) element.numSide = 8;
    else{
        gmsh::logger::write("The element name is not supported by the program.", "error");
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
    // It is the number of edge nodes divided by the number of elements.
    element.numberFrontierNode = element.frontierNode.size()/(element.numSide*element.elementTag.size());

    // Gets the shape functions of the elements at the Gauss points;
    gmsh::model::mesh::getBasisFunctions(element.elementType[0], integrationType, "IsoParametric",\
                                element.gaussPointsParam, element.numComponentShapeFunctions,\
                                element.shapeFunctionsParam);

    // Gets the gradient of the shape functions of the elements at the Gauss points;
    gmsh::model::mesh::getBasisFunctions(element.elementType[0], integrationType, "GradIsoParametric",\
                                            bin1, element.numComponentShapeFunctionsGrad,\
                                            element.shapeFunctionsGradParam);

    // Gets the jacobian at each Gauss point of the elements.
    gmsh::model::mesh::getJacobians(element.elementType[0], integrationType, bin2, element.jacobians,\
                                    bin3);

}