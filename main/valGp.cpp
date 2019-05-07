/* This functions computes the approximated function (Galerkin) at the Gauss points of the elements
   described by "element" in parametric coordinates. It stores the results in the "result" vector in the form
   [e1G1, e1G2, ... , e2G1, e2G2, ...]. */

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void valGp(Quantity & u, const Element & mainElement, const Element & frontierElement, int numU){

    std::size_t i, j, k, l;

    std::fill(u.gp.begin(), u.gp.end(), std::make_pair(0,0));

    for(i = 0; i < frontierElement.elementTag.size(); ++i) // Loop over the elements
        for(j = 0; j < frontierElement.numGp; ++j) // Loop over the Gauss Points
            for(k = 0; k < frontierElement.numNodes; ++k) // loop over the nodes
                for(l = 0; l < numU; ++l)
                {
                    int gpIndex    = i * frontierElement.numGp * numU + j * numU + l;
                    int shapeIndex = j * frontierElement.numNodes + k;
    
                    int frontNodeIndex = i * frontierElement.numNodes + k;

                    int mainNodeIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes * numU + \
                                         frontierElement.nodeCorrespondance[frontNodeIndex].first * numU + l;

                    int mainNodeIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes * numU + \
                                         frontierElement.nodeCorrespondance[frontNodeIndex].second * numU + l;
                    

                    u.gp[gpIndex].first += u.node[mainNodeIndex1] * \
                                           frontierElement.shapeFunctionsParam[shapeIndex];
                    

                    if(frontierElement.neighbours[i].second >= 0)
                        u.gp[gpIndex].second += u.node[mainNodeIndex2] * \
                                                frontierElement.shapeFunctionsParam[shapeIndex];
                    
                    else if(frontierElement.neighbours[i].second < -1)
                        u.gp[gpIndex].second += u.bound[mainNodeIndex1] * \
                                                frontierElement.shapeFunctionsParam[shapeIndex];
                    

                }

            

}