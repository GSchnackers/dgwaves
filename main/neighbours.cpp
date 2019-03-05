// function that associates to an edge the tags of the two elements connected to that edge.

// INPUT: - nodeTags, the tag of the nots by elemnts as given by getElementByType.
//        - nodeNumber, the number of nodes per element.
//        - elementTags, the tags of elements as given bygetElementByType.
//        - nodes, the vector containing the sorted nodes.
//        - neighbourhood, the vector in which the neigbours are placed.
//          It is of the form [edge1tagneighbour1, edge1tagneghbour2, edge2tagneighbour1, ...].
//        - It has twice the size of nodes and must be allocated BEFORE calling the function.
// OUTPUT: none

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

void neighbours(const std::vector<int> nodeTags, const int nodeNumber,\
               const std::vector<int> elementTags, std::vector<int> nodes,\
               std::vector<int> & neighbourhood){

    // The "i" loop stands for each pair of nodes (in fact each edge).
    for(std::size_t i = 0; i < nodes.size(); i += 2){

<<<<<<< HEAD
        std::vector<int> tmp(2,-1); // Temporary vector containing the tags of the neighbouring elements.
                                    // Initialized at -1 (i.e. no neighbours state). 
=======
        std::vector<int> tmp(2,-1); // Temporary vector containing the tags of the neighbouring elements. 
>>>>>>> 28274124b591580bdd9310cf9a19343e1a3dc581
        int l = 0; 

        // The "j" loop stands for each element.
        for(std::size_t j = 0; j < nodeTags.size(); j += nodeNumber){

            int verificator = 0; // Variable which is incremented each time a node is common to two
                                 // elements.

            // The "k" loop stands for each node of one element.
            // At each iteration, 
            for(std::size_t k = 0; k < nodeNumber; ++k)
                if(nodeTags[j + k] == nodes[i] ||nodeTags[j + k] == nodes[i + 1])
                    ++verificator;

            if(verificator == 2){

                tmp[l] = elementTags[j/nodeNumber];
                ++l;

            }  

        }

        neighbourhood[i] = tmp[0];
        neighbourhood[i + 1] = tmp[1];

    }

}