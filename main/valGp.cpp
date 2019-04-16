/* This functions computes the approximated function (Galerkin) at the Gauss points of the elements
   described by "element" in parametric coordinates. It stores the results in the "result" vector in the form
   [e1G1, e1G2, ... , e2G1, e2G2, ...]. */

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void valGp(Quantity & q, const Element & mainElement, const Element & frontierElement, int compo){

    std::size_t i, j, k, l;

    std::fill(q.gp.begin(), q.gp.end(), std::make_pair(0,0));

    for(i = 0; i < frontierElement.elementTag.size(); ++i) // Loop over the elements
        for(j = 0; j < frontierElement.numGp; ++j) // Loop over the Gauss Points
        { 
            int gpIndex = i * frontierElement.numGp + j;
            for(k = 0; k < frontierElement.numNodes; ++k) // loop over the nodes
            {
                
                int shapeIndex = j * frontierElement.numNodes + k;
   
                int frontNodeIndex = i * frontierElement.numNodes + k;

                int mainNodeIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes + \
                                    frontierElement.nodeCorrespondance[frontNodeIndex].first;

                int mainNodeIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes + \
                                    frontierElement.nodeCorrespondance[frontNodeIndex].second;
                

                q.gp[gpIndex].first += q.node[mainNodeIndex1] * \
                                        frontierElement.shapeFunctionsParam[shapeIndex];
                

                if(frontierElement.neighbours[i].second >= 0)
                    q.gp[gpIndex].second += q.node[mainNodeIndex2] * \
                                        frontierElement.shapeFunctionsParam[shapeIndex];
                
                else
                {
                    q.gp[gpIndex].second += q.bound[mainNodeIndex1] * \
                                            frontierElement.shapeFunctionsParam[shapeIndex];
                    //std::cout << q.gp[gpIndex].second << std::endl;
                }
                

            }

            //std::cout << q.gp[gpIndex].first << " " << q.gp[gpIndex].second << std::endl;

        }

}