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

    int entityTag; // stores the entity of the element if it is confined in one particular entity.

    // Those three properties define the true nature of the element.
    std::string name; // Name of the element.
    int dim; // dimension of the element.
    int order; // Order of the element.
    int numSide; // Number of sides frontier elements on one element.

    // Useful properties common to all elements of the element.

    std::vector<double> parametricCoordinates; // Parametric coordinates of the element.

    int numNodes; // Number of nodes of the element.

    std::vector<double> shapeFunctionsParam; // Shape functions at the Gauss points in parametric coordinates.
    int numCompoShape; // Number of shape functions.
    int numCompoShapeGrad; // number of components of the gradient of shape functions.

    std::vector<double> shapeFunctionsGradParam; // Shape functions gradient at the Gauss points in parametric coordinates.
    
    std::vector<double> gaussPointsParam; // Gauss points in parametric coordinates for shape functions.
    std::vector<double> gaussPointsParamGrad; // Gauss points param for gradient of shape functions.
    std::vector<double> gaussPoints; // Gauss points real coordinates.

    int numGp; // Contains the number of gauss points on the element.

    int numberFrontierNode; // The number of nodes per edge of an element.

    std::vector<double> massMatrix; // Mass matrix of the element.
    std::vector<double> massMatrixInverse; // Inverse of the mass matrix of the element.

    std::vector<double> stiffnessMatrixX; // stiffness matrix of the element (x component).
    std::vector<double> stiffnessMatrixY; // stiffness matrix of the element (y component).
    std::vector<double> stiffnessMatrixZ; // stiffness matrix of the element (z component).

    std::vector<double> jacobians; // Jacobians of the element at the gauss points in real coordinates.
    std::vector<double> jacobiansInverse; // Inverse of the real jacobian.
    std::vector<double> jacobiansDet; // Jacobians of the element at the gauss points.

    std::vector<int> frontierNode; // Nodes at the frontier of the element.

    std::vector<int> elementTag; // Tag of all elements.

    std::vector<int> nodeTags; // Tags of the node of each elements. e1N1,e1N2,...,e1Nn, e2N1,...

    std::vector<std::pair<int,int>> neighbours;// Neighbours of the elements. Useful only for the frontier elements.
    std::vector<std::pair<int, int>> nodeCorrespondance; // Vector of pair. nodecorrespondance[i * numNodes + j] contains the indices of the jth neighbour nodes on the ith frontier. firs is for the first element neighbour, second is for the second neighbour.  

    std::vector<double> normals; // Vector containing the normals. Only useful for frontier elements.

};

typedef struct Element Element;

// Structure that deals with the viewing of the results.
struct View{

    std::string name; // The name of the view.
    std::string modelName; // The name of the model attached to view.
    std::string dataType; // The datatype to be put inside the view
    int tag; // The view tag.
    std::vector<std::vector<double>> data; // The data to be viewed.

};

typedef struct View View;

// Deals with the unknowns of the problem and the fluxes.
struct Quantity{

    std::vector<double> node; // Values of the quantity at the nodes of the elements.
    std::vector<std::pair<double, double>> gp; // value of the quantity at the gauss points.
    std::vector<double> direction; // The direction in which the information propagates.
    std::vector<double> bound; // set of 0 and 1 value. It is 1 if the corresponding node is at the boundary, it is 0 if not.
    std::vector<double> num; // numerical value
    std::vector<double> boundSign; // Contains the type of the boundary condition applied to the node.

};

typedef struct Quantity Quantity;

#endif