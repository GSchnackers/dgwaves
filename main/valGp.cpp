/* This functions computes the approximated function (Galerkin) at the Gauss points of the elements
   described by "element" in parametric coordinates. It stores the results in the "result" vector in the form
   [e1G1, e1G2, ... , e2G1, e2G2, ...]. */

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void valGp(const std::vector<double> u, const Element & element, std::vector<double> & result){

    std::size_t i, j, k;

    std::fill(result.begin(), result.end(), 0);

    for(i = 0; i < element.elementTag.size(); ++i) // loop over the elements
        for(j = 0; j < element.numGp; ++j) // loop over the Gauss Points
            for(k = 0; k < element.numNodes ; ++k) // Loope over the nodes (i.e. the shape functions) of the element.
            {

                int resIndex = i * element.numGp + j; // index of the result.
                int uIndex = i * element.numNodes + k; // index of the u vector.
                int shapeIndex = k * element.numGp + j; // index of the gauss point coordinates.

                result[resIndex] += u[uIndex] * element.shapeFunctionsParam[shapeIndex];

            }

}