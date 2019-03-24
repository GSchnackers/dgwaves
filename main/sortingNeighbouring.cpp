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
    std::size_t i, j, k;

    // Temporary variables.
    std::vector<std::pair<int,int>> tmpNeighbours(mainElement.numSide * mainElement.elementTag.size());
    std::vector<int> tmpSorted(tmpNeighbours.size() * mainElement.numberFrontierNode);

    // Counter
    int countNeigh = 0, countNode = 0;

    // Useful variables.
    int totalNumberFrontierNode = mainElement.numSide * mainElement.numberFrontierNode;

    // Initialisation of the node sorted vector.
    for(i = 0; i < mainElement.numberFrontierNode; ++i){

        tmpSorted[i] = mainElement.frontierNode[i];
        ++countNode;
    }

    // Neighbours initializations.
    tmpNeighbours[0].first = 0;
    tmpNeighbours[1].first = 1;
    ++countNeigh;

    // Loop over the frontier nodes. i is the index of the frontier nodes.
    for (i = mainElement.numberFrontierNode; i < mainElement.frontierNode.size(); i += mainElement.numberFrontierNode)
    {
        int elemIndexi = i/totalNumberFrontierNode; // Index of an element with respect to i
        int gradIndexi = 3*elemIndexi; // Index of the gradient.
        int edgeEndi = i + mainElement.numberFrontierNode - 1; // index of the end of the edge (in the unsorted vector).

        // Loop over the sorted nodes.
        for (j = 0; j < tmpSorted.size(); j += mainElement.numberFrontierNode)
        {
            int frontierIndexj = j/mainElement.numberFrontierNode; // Index of a frontier element with respect to j.
            int edgeEndj = j + mainElement.numberFrontierNode - 1; // index of the end of the edge in the sorted vector.

            bool condition = (mainElement.frontierNode[i] == tmpSorted[edgeEndj] \
                            && mainElement.frontierNode[edgeEndi] == tmpSorted[j]) || \
                            (mainElement.frontierNode[i] == tmpSorted[j] \
                            && mainElement.frontierNode[edgeEndi] == tmpSorted[edgeEndj]); // check if two edges coincide. Two edges coincide only in opposition, that is the first node of the first edge equals the last node of the second edge etc.


            // If the edge is already in the neighbour, the element INDEX (not tag) is then added
            // to the neighbourhood of the edge.
            if(condition)
            {
                
                if(tmpNeighbours[frontierIndexj].first != elemIndexi)
                    tmpNeighbours[frontierIndexj].second = elemIndexi;
                
                break;
            }

            // If the end of the list is reached without finding a match to the edge in the 
            // nodeSorted vector, the edge is added to the nodeSorted vector. A neighbour
            // entry is created at the same time and filled with the index of the element
            // and -1, assuming the element has only one neighbour at first. The normal to
            // the edge is also created at the same time.

            if(j + mainElement.numberFrontierNode == tmpSorted.size())
            {
                for(k = 0; k < mainElement.numberFrontierNode; ++k){

                    tmpSorted[countNode] = mainElement.frontierNode[i + k];
                    ++countNode;

                } 

                tmpNeighbours[countNeigh].first = elemIndexi;
                tmpNeighbours[countNeigh].second = -1;
                ++countNeigh;
                
            }
        }
    }

    frontierElement.neighbours.resize(countNeigh);
    nodeSorted.resize(countNode);
    
    for(i = 0; i < countNeigh; ++i) frontierElement.neighbours[i] = tmpNeighbours[i];
    for(i = 0; i < countNode; ++i) nodeSorted[i] = tmpSorted[i];

}