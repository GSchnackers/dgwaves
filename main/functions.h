#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "structures.h"


// Initialization of the properties of the element of a certain dim and a certain type with a number of Gauss points given by GaussType.
void Initialization(Element & element, const int meshDim, std::string integrationType, bool frontier = false);

// Creation of the frontier of the main elements. Each frontier is only counted once.
void frontierCreation(const Element mainElement, Element & frontierElement, const int meshDim,\
                      const std::vector<int> sortedNodes);

// Sorting the nodes of the list of the edges to remove duplicate and giving the neighbours of each edge.
void sortingNeighbouring(const Element & mainElement, Element & frontierElement,\
                         std::vector<int> & nodeSorted);

// Functions that allows to inverse a matrix.
void invert(std::vector<double> matrix, std::vector<double> & inverse);

// Function that inverses the jacobian.
void getJacobiansInverse(Element & element);

// Gets the gradient in real coordinates of the shape function of each elements represented by element.
void getRealGradient(Element & element);

// Gets the normal to all edge elements.
void normals(Element & frontierElement);

// Functions that computes the main matrices of the DG fem, that is the mass and the stiffiness matrices.
void matrixMaker(Element & element, std::string matrixType);


/*// Compute the integration of Gauss
void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & functions,
 const std::vector<double> & determinants, std::vector<double> & matrix,
  const int numElements, const int numGaussPoints, const int numNodes); */

#endif