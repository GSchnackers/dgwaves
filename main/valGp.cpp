/* This functions computes the approximated function (Galerkin) at the Gauss points of the elements
   described by "element" in parametric coordinates. It stores the results in the "result" vector in the form
   [e1G1, e1G2, ... , e2G1, e2G2, ...]. */

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void valGp(Quantity & u, const Element & mainElement, const Element & frontierElement){

    std::size_t i, j, k;

    u.gp.resize(frontierElement.elementTag.size() * frontierElement.numGp);

    std::fill(u.gp.begin(), u.gp.end(), 0);

    for(i = 0; i < frontierElement.elementTag.size(); ++i) // Loop over the elements
        for(j = 0; j < frontierElement.numGp; ++j) // Loop over the Gauss Points
            for(k = 0; k < frontierElement.numNodes ; ++k) // Loop over the nodes (i.e. the shape functions) of the element.
            {

                int indexGp = i * frontierElement.numGp + j;
                int indexNode1 = frontierElement.neighbours[i].first * frontierElement.numNodes + k;
                int indexShape = j * frontierElement.numNodes + k;

                u.gp[indexGp].first += u.node[indexNode1] * frontierElement.shapeFunctionsParam[indexShape];

                // Case we do not have a frontier.
                if(frontierElement.neighbours[i].second > -1)
                {
                    int indexNode2 = frontierElement.neighbours[i].second * frontierElement.numNodes + k;
                    
                    u.gp[indexGp].second += u.node[indexNode2] * frontierElement.shapeFunctionsParam[indexShape];

                }

                // In case we have a frontier.
                else if(frontierElement.neighbours[i].second == -1)
                    u.gp[indexGp].second += u.node[indexNode1] * frontierElement.shapeFunctionsParam[indexShape];

            }

}