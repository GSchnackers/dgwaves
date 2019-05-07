/*

This file contains the fuction that stores the indices corresponding to two neighbour nodes of the frontier 
elements.

*/

#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include "structures.hpp"


// This function links the nodes of the frontier elements with their indices in the global numerotation.
void correspondance(const Element & mainElement, Element & frontierElement){

    std::size_t i, j, k;
    frontierElement.nodeCorrespondance.resize(frontierElement.elementTag.size() * frontierElement.numNodes);

    for(i = 0; i < frontierElement.elementTag.size(); ++i) // loop over the frontier elements
        for(j = 0; j < frontierElement.numNodes; ++j) // loop over the nodes of one frontier element.
        {
            int frontNodeIndex = i * frontierElement.numNodes + j; // index of the correspondance vector case to be filled.

            for(k = 0; k < mainElement.numNodes; ++k) // Seraching for the first index of the first neighbour
            {
                
                int neighIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes + k;
                
                if(frontierElement.nodeTags[frontNodeIndex] == mainElement.nodeTags[neighIndex1])
                {
                    frontierElement.nodeCorrespondance[frontNodeIndex].first = k;
                    break;
                }

            }

            if(frontierElement.neighbours[i].second > 0) // In case there is a second neighbour, searches the index of the second neighbour node.
                for(k = 0; k < mainElement.numNodes; ++k)
                {
                    int neighIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes + k;

                    if(frontierElement.nodeTags[frontNodeIndex] == mainElement.nodeTags[neighIndex2])
                    {
                        frontierElement.nodeCorrespondance[frontNodeIndex].second = k;
                        break;
                    }
                    
                }

            else // In case there is no neighbour, simply fills with -1.
                frontierElement.nodeCorrespondance[frontNodeIndex].second = -1;
            
        }
}