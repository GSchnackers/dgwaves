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
            int gpIndex = i * frontierElement.numGp + j;
            u.numGp[gpIndex].first = u.numGp[gpIndex].second = 0;
            
            for(k = 0; k < frontierElement.numNodes; ++k)
            { 
                
                int frontNodeIndex = i * frontierElement.numNodes + k;
                int mainNodeIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes + \
                                     frontierElement.nodeCorrespondance[frontNodeIndex].first;
                int mainNodeIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes + \
                                     frontierElement.nodeCorrespondance[frontNodeIndex].second;
                int shapeIndex = j * frontierElement.numNodes + k;
                

                u.numGp[gpIndex].first += u.node[mainNodeIndex1] * \
                                          frontierElement.shapeFunctionsParam[shapeIndex];

                if(frontierElement.neighbours[i].second < 0)
                    u.numGp[gpIndex].second = u.numGp[gpIndex].first;

                else
                    u.numGp[gpIndex].second += u.node[mainNodeIndex2] * \
                                          frontierElement.shapeFunctionsParam[shapeIndex];
                
            }


        }

}