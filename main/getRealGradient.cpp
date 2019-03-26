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

    // Resize of the shap function gradient vector.
    element.shapeFunctionsGrad.resize(element.shapeFunctionsGradParam.size() * element.elementTag.size());
    
    // Multiplication of the isoparametric shape functions gradient by the inverse of the jacobian.
    for(i = 0; i < element.elementTag.size(); ++i) // Loop over the elements. 
        for(j = 0; j < element.numGp; ++j) // Loop on the nodes.
            for(k = 0; k < element.numNodes; ++k) // Loop over the Gauss Points.
            {
                double tmp = 0;
                int m = 0;
                for(l = 0; l < 9; ++l) // Loop over the jacobian.
                {
                    int indexGradParam = j * element.numNodes * element.numCompoShapeGrad +\
                        k * element.numCompoShapeGrad + l % 3;
                    int indexJacob = i * element.numGp * 9 + j * 9 + l;

                    int indexGrad = i * element.numGp * element.numNodes * element.numCompoShapeGrad +\
                                    j * element.numNodes * element.numCompoShapeGrad + \
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