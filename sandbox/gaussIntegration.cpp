#include <cstdio>
#include <iostream>
#include <vector>
#include <gmsh.h>
#include "functions.h"

void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & basisFunctions,
 const std::vector<double> & determinants, std::vector<double> & massMatrix, const int numElements, const int numGaussPoints)
{
    //numGaussPoints = integrationPoints.size()/numElements;
    for(std::size_t i = 0; i < numElements; i++)
    {
        for(std::size_t j = 0; j < 3; j++)
        {
            for(std::size_t k = 0; k < 3; k++)
            {
                double gaussSum = 0;

                for(std::size_t l = 0; l < numGaussPoints; l++)
                {
                    gaussSum += integrationPoints[3 + 4*l] * basisFunctions[3*l + k] * basisFunctions[3*l + j]
                                * determinants[numGaussPoints*i + l];
                }
                    massMatrix[9*i + 3*j + k] = gaussSum;
            }
        }
    }
}