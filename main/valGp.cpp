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

    for(i = 0; i < frontierElement.elementTag.size(); ++i) // Loop over the elements
        for(j = 0; j < frontierElement.numGp; ++j) // Loop over the Gauss Points
        { 
            int frontGpIndex = i * frontierElement.numGp + j; // index of the gauss points at the frontier elements.
            u.numGp[frontGpIndex] = std::make_pair(0,0);

            for(k = 0; k < frontierElement.numNodes ; ++k) // Loop over the nodes (i.e. the shape functions) of the element.
            {

                
                int frontNodeIndex = i * frontierElement.numNodes + k; // index of the nodes on the frontier elements.

                int mainNodes1 = frontierElement.neighbours[i].first * mainElement.numNodes + \
                                 frontierElement.nodeCorrespondance[frontNodeIndex].first; 

                // First neighbour.
                u.numGp[frontGpIndex].first += u.node[mainNodes1] * mainElement.shapeFunctionsParam[j];

                // If there is another neighbour to the considered frontier element.
                if(frontierElement.neighbours[i].second >= 0)
                {
                    // Index of the second node in general numbering.
                    int mainNodes2 = frontierElement.neighbours[i].second * mainElement.numNodes + \
                                    frontierElement.nodeCorrespondance[frontNodeIndex].second; 

                    u.numGp[frontGpIndex].second += u.node[mainNodes2] * mainElement.shapeFunctionsParam[j];
                }
                
                // Case we encounter a boundary condition.
                else
                {
                    u.numGp[frontGpIndex].second += u.bound[frontNodeIndex] * \
                                                    mainElement.shapeFunctionsParam[j];
                }
                

            }
        }

}