// INPUT: -  integrationPoints: contains the parametric coordinates and the weight for each integration point
//        -  functions: vector containing the values of the function to integrate at the gauss points. 
//        -  determinants: contains the determinant of the Jacobian matrix at each Gauss point.
//        -  numElements: number of elements in the mesh.
//        -  numGaussPoints: number of Gauss points in element.
//        -  numNodes: number of nodes of each element.
// OUTPUT: - matrix: the matrix of integrated functions.

#include <cstdio>
#include <iostream>
#include <vector>
#include <gmsh.h>
#include "functions.h"

void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & functions,
 const std::vector<double> & determinants, std::vector<double> & matrix,
  const int numElements, const int numGaussPoints, const int numNodes)
{
    //numGaussPoints = integrationPoints.size()/numElements;
    for(std::size_t e = 0; e < numElements; e++)
        for(std::size_t i = 0; i < numNodes; i++)
            for(std::size_t j = 0; j < numNodes; j++)
            {
                double gaussSum = 0;

                for(std::size_t g = 0; g < numGaussPoints; g++)
                    // function = [e1i1j1g1, e1i1j1g2, ..., e1i1j1gG, e1i1j2g1, ..., e1i2j1g1, ..., e2i1j1g1, ...]
                    // functions[numNodes*numNodes*numGaussPoints*e + numNodes*numGaussPoints*i + numGaussPoints*j + g]
                    gaussSum += functions[numGaussPoints*(numNodes*(numNodes*e + i) + j) + g]
                                * integrationPoints[3 + 4*g] * determinants[numGaussPoints*e + g];
                matrix.push_back(gaussSum);
            }
}