#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include <Eigen/Dense>

void MatrixInverseM(const std::vector<double> shapeFunctions, const std::vector<double> weight,\
                    const int numNodes, const std::vector<double> jacobian,\
                    const int integrationtype, std::vector<double> & massMatrix){
  
    std::size_t i, j, k;
    Eigen::Matrix<double, numNodes, numNodes> tmp = Eigen::Matrix<double, numNodes, numNodes>::Zero();

    // The mass matrix is symmetric.

    for(i = 0; i < numNodes * integrationType; i += integrationType)
        for(j = i; j < numNodes * integrationType; j += integrationType)
            for(k = 0; k < integrationType; ++k){

                tmp(i/integrationType,j/integrationType) += weight[k] * shapeFunction[i + k]\
                                                            * shapeFunction[j + k] * jacobian[k];
                tmp(j/integrationType,i/integrationType) =  tmp(i/integrationType,j/integrationType);

            }
    
    tmp = tmp.inverse();

    for(i = 0; i < numNodes; ++i)
        for(j = i; j < numNodes; ++j)
            for(k = 0; k < integrationType; ++k) massMatrix[numNodes * i + j] = \
                                                 massMatrix[numNodes * j + i] = tmp(i,j);                                                     

}