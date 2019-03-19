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

    // Useful variables.
    int totalNumberFrontierNode = mainElement.numSide*mainElement.numberFrontierNode;

    // Initialisation of the node sorted vector.
    for(i = 0; i < mainElement.numberFrontierNode; ++i)
        nodeSorted.push_back(mainElement.frontierNode[i]);

    std::pair<int, int> pairTmp(0, -1);
    frontierElement.neighbours.push_back(pairTmp);

    // Loop over the frontier nodes.
    for (i = mainElement.numberFrontierNode; i < mainElement.frontierNode.size(); i += mainElement.numberFrontierNode)

        // Loop over the sorted nodes.
        for (j = 0; j < nodeSorted.size(); j += mainElement.numberFrontierNode)
        {
            // Check if both edges.
            bool condition = (mainElement.frontierNode[i] == nodeSorted[j] \
                             && mainElement.frontierNode[i + mainElement.numberFrontierNode - 1]\
                              == nodeSorted[j + mainElement.numberFrontierNode - 1])\
                            || (mainElement.frontierNode[i]\
                             == nodeSorted[j + mainElement.numberFrontierNode - 1] \
                            && mainElement.frontierNode[i + mainElement.numberFrontierNode - 1]\
                             == nodeSorted[j]);

            // If the edge is already in the neighbour, the element INDEX (not tag) is then added
            // to the neighbourhood of the edge.
            if(condition){

                if(frontierElement.neighbours[j/(mainElement.numberFrontierNode)].first\
                    != i/totalNumberFrontierNode)
                    frontierElement.neighbours[j/(mainElement.numberFrontierNode)].second = \
                                        i/totalNumberFrontierNode;
                break;

            }

            // If the end of the list is reached without finding a match to the edge in the 
            // nodeSorted vector, the edge is added to the nodeSorted vector. A neighbour
            // entry is created at the same time and filled withe the index of the element
            // and -1, assuming the element has only one neighbour at first. A normal entry is also
            // Created at the same time.

            if(j + mainElement.numberFrontierNode == nodeSorted.size())
            {
                for(k = 0; k < mainElement.numberFrontierNode; ++k) 
                    nodeSorted.push_back(mainElement.frontierNode[i + k]);

                std::pair<int, int> pairTmp(i/totalNumberFrontierNode, -1);

                frontierElement.neighbours.push_back(pairTmp);

                if(mainElement.dim == 2)
                    for(k = 0; k < frontierElement.shapeFunctionsGradParam.size(); k += 3){

                        double norm = sqrt(frontierElement.shapeFunctionsGradParam[k + 1] *\
                                           frontierElement.shapeFunctionsGradParam[k + 1] +\
                                           frontierElement.shapeFunctionsGradParam[k] *
                                           frontierElement.shapeFunctionsGradParam[k]);

                        frontierElement.normals.push_back(\
                                                frontierElement.shapeFunctionsGradParam[k + 1]/norm);
                        frontierElement.normals.push_back(\
                                                -frontierElement.shapeFunctionsGradParam[k]/norm);
                        frontierElement.normals.push_back(0.);

                    }
                else if(mainElement.dim == 3){

                    for(k = 0; k < frontierElement.elementTag.size(); ++k){

                        std::vector<double> normalTmp;
                        std::vector<double> gaussTmp = {frontierElement.gaussPointsParam[3*k],\
                                                        frontierElement.gaussPointsParam[3*k + 1]};
                        gmsh::model::getNormal(frontierElement.elementTag[k], gaussTmp, normalTmp);
                        for(std::size_t l = 0; l < normalTmp.size(); ++l)
                            frontierElement.normals.push_back(normalTmp[l]);

                    }

                }

                
            }

        }

}