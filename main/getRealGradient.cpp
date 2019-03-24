#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <Eigen/Dense>
#include <cmath>
#include "functions.h"

// Functions that computes the real gradient of the elements on the basis of the parametric gradients.

void getRealGradient(Element & element)
{
    std::size_t i, j, k, l;

    // Temporary vector for shape functions.
    
    element.shapeFunctionsGrad.resize(element.shapeFunctionsGradParam.size() * element.elementTag.size());
    std::cout << element.shapeFunctionsGrad.size() << std::endl;
    
    // Multiplication of the isoparametric shape functions gradient by the inverse of the jacobian.
    for(i = 0; i < element.elementTag.size(); ++i) // Loop over the elements. 
        for(j = 0; j < element.numNodes; ++j) // Loop on the nodes.
            for(k = 0; k < element.numGp; ++k) // Loop over the Gauss Points.
            {
                double tmp = 0;
                int m = 0;
                for(l = 0; l < 9; ++l) // Loop over the jacobian.
                {
                    int indexGradParam = j * element.numGp * element.numCompoShapeGrad +\
                        k * element.numCompoShapeGrad + l % 3;
                    int indexJacob = i * element.numGp * 9 + k * 9 + l;

                    int indexGrad = i * element.numNodes * element.numGp * element.numCompoShapeGrad +\
                                    j * element.numGp * element.numCompoShapeGrad + \
                                    k * element.numCompoShapeGrad;

                    tmp += element.jacobiansInverse[indexJacob] * element.shapeFunctionsGradParam[indexGradParam];

                    if(l % 3 == 2)
                    {
                        element.shapeFunctionsGrad[indexGrad + m] = tmp;
                        ++m;
                        tmp = 0;
                    }
                }
            } 

}