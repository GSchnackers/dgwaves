// INPUT: -  integrationPoints, the points where the gauss integration is performed + the weights of each point.
//        -  functions: vector containing the values of the function to integrate at the gauss points. 
//        -  determinants: contains the jacobian at each gauss point.
//        -  numElements: The number of elements in the mesh.
//        -  numGaussPoints: the number of Gauss points in an element.
// OUTPUT: - matrix: the matrix of integrated functions.

#include <cstdio>
#include <iostream>
#include <vector>
#include <gmsh.h>
#include "functions.h"

void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & functions,
 const std::vector<double> & determinants, std::vector<double> & matrix, const int numElements, const int numGaussPoints)
{
    //numGaussPoints = integrationPoints.size()/numElements;
    for(std::size_t i = 0; i < numElements; i++)
        for(std::size_t j = 0; j < 3; j++)
            for(std::size_t k = 0; k < 3; k++)
            {
                double gaussSum = 0;

                for(std::size_t l = 0; l < numGaussPoints; l++)
                    gaussSum += integrationPoints[3 + 4*l] * functions[3*l + k] * determinants[numGaussPoints*i + l];
                
                matrix[9*i + 3*j + k] = gaussSum;
            }
        
    
}