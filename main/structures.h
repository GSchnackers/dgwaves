#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"


// Structure Element represent the common characteristic of element of dimension dim.

struct Element{

        std::string name; // Name of the element.
        std::vector<double> parametricCoordinates; // Parametric coordinates of the element.
        std::vector<int> elementTypes; // Type of element.
        int dim; // dimension of the element.
        int order; // Order of the element.
        int numNodes; // Number of nodes of the element.

};

struct Entity{

    int entityTag;
    std::vector<int> nodeTags;
    std::vector<int> elementEdgeNodes;
    std::vector<int> edges;
    std::vector<int> elementTag;
    std::vector<int> neighbours;

};
#endif