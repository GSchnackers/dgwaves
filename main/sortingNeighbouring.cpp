#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include "functions.h"
#include "structures.h"

// This functions takes the nodes on the frontier of the main elements and returns the frontier
// element nodes with no doubles in the form [e1n1, e2n2, ..., e1nN, e2n1,...] and the neighbours
// in form of a vector pair [e1 (neighbour1Index, neighbour2Index),
// e2 (neighbour1Index, neighbour2Index), ...]. The first element always correspond to an element
// neighbouring the frontier, while the other is -1 if there is no other neighbour to that element.

void sortingNeighbouring(const Element & mainElement, Element & frontierElement,\
                         std::vector<int> & nodeSorted)
{
    // Increment variables.
    std::size_t i, j, k, l;

    // Temporary variables.
    std::vector<std::pair<int,int>> tmpNeighbours(mainElement.numSide * mainElement.elementTag.size());
    std::vector<int> tmpSorted(tmpNeighbours.size() * mainElement.numberFrontierNode);

    // Counter
    int sizeNeigh = 0, sizeNode = 0;

    // Useful variables.
    int totalNumberFrontierNode = mainElement.numSide * mainElement.numberFrontierNode;

    // Initialisation of the node sorted vector.
    for(i = 0; i < mainElement.numberFrontierNode; ++i){

        tmpSorted[i] = mainElement.frontierNode[i];
    }

    sizeNode += mainElement.numberFrontierNode;

    // Neighbours initializations.
    tmpNeighbours[0].first = 0;
    tmpNeighbours[0].second = -1;
    ++sizeNeigh;

    // Loop over the frontier nodes. i is the index of the frontier nodes.
    for(i = mainElement.numberFrontierNode; i < mainElement.frontierNode.size(); i += mainElement.numberFrontierNode)
    { 
        int elementIndexi = i/(mainElement.numSide * mainElement.numberFrontierNode);
        for (j = 0; j < sizeNode; j += mainElement.numberFrontierNode) // Loop over the sorted nodes.
        {
            int countNode = 0;
            int neighbourIndexj = j/mainElement.numberFrontierNode;
            // Comparison of the nodes of the unsorted vector and the sorted vector. If each node of the
            // unsorted find a match with each node of the unsorted, then the frontier has already been
            // taken into accound and is ignored. The second neighbour is set at that stage.
            // If it is not the case, then the new frontier is added to the sorted vector.

            for(k = 0; k < mainElement.numberFrontierNode; ++k)
                for(l = 0; l < mainElement.numberFrontierNode; ++l)
                    if(mainElement.frontierNode[i + k] == tmpSorted[j + l])
                    {
                        ++countNode;
                        break;
                    }

            if(countNode == mainElement.numberFrontierNode)
            {   
                if(tmpNeighbours[neighbourIndexj].second != elementIndexi)  
                    tmpNeighbours[neighbourIndexj].second = elementIndexi;
                break;
            }

            else if(j == sizeNode - mainElement.numberFrontierNode)
            {
                
                for(k = 0; k < mainElement.numberFrontierNode; ++k)
                {
                    
                    tmpSorted[sizeNode] = mainElement.frontierNode[i + k];
                    ++sizeNode;
                }

                tmpNeighbours[sizeNeigh].first = elementIndexi;
                tmpNeighbours[sizeNeigh].second = -1;
                
                ++sizeNeigh;
                break;
            }
        }
    }

    frontierElement.neighbours.resize(sizeNeigh);
    nodeSorted.resize(sizeNode);
    

    // Final definiton of the neighbour vector.
    for(i = 0; i < sizeNeigh; ++i)
        frontierElement.neighbours[i] = tmpNeighbours[i];

    // Final definition of the sorted node vector.
    for(i = 0; i < sizeNode; ++i)
        nodeSorted[i] = tmpSorted[i];

}