#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

// Structure Element represent the common characteristic of element of dimension dim.

struct Element{

    // Type of element, its true nature
    std::vector<int> elementType;

    // Those three properties define the true nature of the element.
    std::string name; // Name of the element.
    int dim; // dimension of the element.
    int order; // Order of the element.

    // Useful properties common to all elements of the element.

    std::vector<double> parametricCoordinates; // Parametric coordinates of the element.

    int numNodes; // Number of nodes of the element.

    std::vector<double> shapeFunctionsParam; // Shape functions at the Gauss points in parametric coordinates.
    int numComponentShapeFunctions; // Number of components of each shapeFunctions.

    std::vector<double> shapeFunctionsGradParam; // Shape functions gradient at the Gauss points in parametric coordinates.
    int numComponentShapeFunctionsGrad; // Number of components of each gradient of shapeFunctions.
    
    std::vector<double> gaussPointsParam; // Gauss points in parametric coordiantes.

    int numberEdgeNode; // The number of nodes per edge of an element.

    std::vector<double> massMatrix; // Mass matrix of the element.
    std::vector<double> stiffnessMatrix; // stiffness matrix of the element.

    std::vector<double> jacobians; // Jacobians of the element at the gauss points.

    std::vector<int> frontierNode; // Nodes at the frontier of the element.

    std::vector<int> elementTag; // Tag of all elements.

    std::vector<int> nodeTags; // Tags of the node of each elements. e1N1,e1N2,...,e1Nn, e2N1,...

};

typedef struct Element Element;

#endif