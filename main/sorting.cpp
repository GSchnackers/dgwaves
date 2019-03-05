#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"

void sorting(std::vector<int> & nodes)
{
    std::vector<int> sortingNodes;
    sortingNodes.push_back(nodes[0]);
    sortingNodes.push_back(nodes[1]);
    for (std::size_t i = 0; i < nodes.size(); i += 2)
    {
        for (std::size_t j = 0; j < sortingNodes.size(); j += 2)
        {
            // Check if edges is already in sortingNodes in the same direction
            if(nodes[i] == sortingNodes[j] && nodes[i+1] == sortingNodes[j+1])
            {
                j = sortingNodes.size();
            }
            // Check if edges is already in sortingNodes in the opposite direction
            else if(nodes[i] == sortingNodes[j+1]  && nodes[i+1] == sortingNodes[j])
            {
                j = sortingNodes.size();
            }

            // If the edge is not already in sortingNodes, we add it
            if(j+2 == sortingNodes.size())
            {
                sortingNodes.push_back(nodes[i]);
                sortingNodes.push_back(nodes[i+1]);
            }
        }
    }

    nodes = sortingNodes;
}