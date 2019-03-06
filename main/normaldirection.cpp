// Functions that returns a vector containing binary value. It gives info about the direction
// (outward or inward) of the normal with respect to the element.

// INPUT: -nodes, the sorted nodes.
//        -normal, the normal to the edges defines by the pair of nodes in nodes.
//        -neighbours, the vector containing the neighbours sharing the same edge.
//        -nodeCoords, the coordinates of the nodes of the elements.
//        -outward, the vector that will contain the binary value. If it is -1, it is inward
//        if it is 0, there is no element on the other side of the edge to point to, 1 it is
//        an outward pointing normal.

// OUTPUT: NONE.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

void normaldirection(const std::vector<int> nodes, const std::vector<double> normal,\
                     const std::vector<int> neighbours, const std::vector<double> nodeCoords,\
                     std::vector<int> & outward){

    std::vector<int> initialize(nodes.size(), 0);
    std::vector<double> tmp(2,0);
    outward = initialize;

    // Run through each pair of neighbours of each edge.
    for(std::size_t i = 0; i < neighbours.size(); i += 2)

        // Also run through the neighbours to find the edges that have an element in common.
        for(std::size_t j = 0; j < neighbours.size(); j += 2){

            int k = -1, l = -1;
            // Check all possible case of compatibility. 
            if(neighbours[i] == neighbours[j] && nodes[i] != nodes[j] && neighbours[i] != -1){

                k = i; l = j; 
            
            }

            else if(neighbours[i] == neighbours[j + 1] && nodes[i] != nodes[j + 1] && neighbours[i] != -1){

                k = i; l = j + 1;

            }

            else if(neighbours[i + 1] == neighbours[j] && nodes[i + 1] != nodes[j] && neighbours[i + 1] != -1){

                k = i + 1; l = j;

            }

            else if(neighbours[i + 1] == neighbours[j + 1] && nodes[i + 1] != nodes[j + 1] && neighbours[i + 1] != -1){

                k = i + 1; l = j + 1;

            }

            if(k && l){

                tmp[0] = -(nodeCoords[3*l+1] - nodeCoords[3*k+1]);
                tmp[1] = nodeCoords[3*l] - nodeCoords[3*k];

                if(tmp[0] * normal[i] + tmp[1] * normal[i + 1] < 0){

                    outward[i] = 1;
                    if(neighbours[i + 1] != -1) outward[i + 1] = -1;

                } 
                else if(tmp[0] * normal[i] + tmp[1] * normal[i + 1] > 0){

                    outward[i] = -1;
                    outward[i + 1] = 1;

                } 

                break;

            }
        
        }

        

    

    

}