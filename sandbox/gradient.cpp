// Computes the divergence of a lagrangian shape function at point x,y,z, no matter the parametric
// or real system of coordinate.

// Input: - f, the value of the lagrangian function (parametric or cartesian coordinates).
//        - pointPosition, contains the real or parametric coordinates of the point where f 
//        is evaluated.
//        - nodePosition contains the real or parametric coordinates of the nodes of the element.

// Output: - grad, the gradient of f at the specified position.

#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"

void gradient(const double f,const  std::vector<double> pointsPositions,\
                             const std::vector<double> nodesPositions, std::vector<double> & grad){

    for(std::size_t i = 0; i < nodesPositions.size(); i += 3)
        for(std::size_t j = 0; j < grad.size() \
            && pointsPositions[i + j] != nodesPositions[i + j] ; ++i)
                grad[i] += f/(pointsPositions[i + j] - nodesPositions[i + j]);
                

}