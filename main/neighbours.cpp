// function that associates to an edge formed by a pair of adjacent nodes the tags of the two elements
// by to that edge.

// INPUT: nodeTags, the tag of the nots by elemnts as given by getElementByType.
//        nodeNumber, the number of nodes per element.
//        elementTags, the tags of elements as given bygetElementByType.
// OUTPUT: nodes, the vector containing the edges and the associated neighbouring elements in the
//         form [edge1n1, edge1n2, edge1tag1, edge1tag2, edge2n1, ...]. If there is only one 
//         neighbour to edgex, the value of edgextag2 = -1.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

void neighbours(const std::vector<int> nodeTags, const int nodeNumber,\
               const std::vector<int> elementTags, std::vector<int> & nodes){

    std::vector<int> neighbourhood(2*nodes.size());

    // The "i" loop stands for each pair of nodes (in fact each edge).
    for(std::size_t i = 0; i < nodes.size(); i += 2){

        std::vector<int> tmp(2,-1); // Temporary vector containing the tags of the neighbouring elements. 
        int l = 0; 

        // The "j" loop stands for each element.
        for(std::size_t j = 0; j < nodeTags.size(); j += nodeNumber){

            int verificator = 0;

            // The "k" loop stands for each node of one element. At each iteration, 
            for(std::size_t k = 0; k < nodeNumber; ++k)
                if(nodeTags[j + k] == nodes[i] ||nodeTags[j + k] == nodes[i + 1])
                    ++verificator;

            if(verificator == 2){

                tmp[l] = elementTags[j/nodeNumber];
                ++l;

            }  

        }

        int m = 2*i;

        neighbourhood[m] = nodes[i];
        neighbourhood[m + 1] = nodes[i + 1];
        neighbourhood[m + 2] = tmp[0];
        neighbourhood[m + 3] = tmp[1];

    }

    nodes = neighbourhood;
}