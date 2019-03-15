#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"


// Structure Element represent the common characteristic of element of dimension dim.

struct ElementProperties{

    // Type of element, its true nature
    std::vector<int> elementType;

    // Those three properties define the true nature of the element.
    std::string name; // Name of the element.
    int dim; // dimension of the element.
    int order; // Order of the element.

    // Useful properties common to all elements of the element.

    std::vector<double> parametricCoordinates; // Parametric coordinates of the element.
    int numNodes; // Number of nodes of the element.
    std::vector<double> shapeFunctionsParam; // Shape functions at the Gauss points in parametric coordiantes.
    std::vector<double> gaussPointsParam; // Gauss points in parametric coordiantes.
    std::vector<double> massMatrix; // Mass matrix of the element.
    std::vector<double> stiffnessMatrix; // stiffness matrix of the element.

};


struct Entity{

    // Tag of the entity.
    int entityTag;

    // Tag of the nodes in the entity.
    std::vector<int> nodeTags;

    struct ElementProperties entityElementProperties1D; // Properties of all 1D elements.
    struct ElementProperties entityElementProperties2D; // Properties of all 2D elements.
    struct ElementProperties entityElementProperties3D; // Properties of all 3D elements.

    // Tags of all elements.
    std::vector<int> elementTag1D;
    std::vector<int> elementTag2D;
    std::vector<int> elementTag3D;

    // Edge nodes of the 2D and 3D elements.
    std::vector<int> elementEdgeNode2D;
    std::vector<int> elementEdgeNode3D;

    // Face nodes of the 3D elements
    std::vector<int> elementFaceNodes3D;

    // Sorted edges nodes.
    std::vector<int> elementEdgeNode2DSorted;
    std::vector<int> elementEdgeNode3DSorted;

    // Face nodes of the 3D elements.
    std::vector<int> elementFaceNodes3DSorted;
    
    // Vector containing the neighbours of the elements.
    std::vector<int> neighbours2D;
    std::vector<int> neighbours3D;

};
#endif