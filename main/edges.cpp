#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"
#include "structures.h"

void edges(const int numNodes, const int numEdgeNode, const std::vector<int> elementEdgeNode,\
           const std::vector<int> elementTag, std::vector<int> & elementEdgeNodeSorted,\
           gmsh::vectorpair & neighbours)
{
    std::size_t i, j, k;

    elementEdgeNodeSorted.push_back(elementEdgeNode[0]);
    elementEdgeNodeSorted.push_back(elementEdgeNode[1]);

    for (i = 0; i < elementEdgeNode.size(); i += numEdgeNode)

        for (j = 0; j < elementEdgeNode.size(); j += numEdgeNode)
        {
            // Check if edges is already in entity.edges in the same or opposite direction.
            bool condition = (elementEdgeNode[i] == elementEdgeNodeSorted[j] \
                             && elementEdgeNode[i + numEdgeNode] == elementEdgeNodeSorted[j + numEdgeNode])\
                            || (elementEdgeNode[i] == elementEdgeNodeSorted[j + numEdgeNode] \
                            && elementEdgeNode[i + numEdgeNode] == elementEdgeNodeSorted[j]);

            if(condition){

                neighbours[j/(numEdgeNode*numNodes)].second = elementTag[j/(numEdgeNode*numNodes)];
                break;

            }

            // If the edge is not already in entity.edges, we add it to the sorting node.
            if(j + numNodes == elementEdgeNode.size())
            {
                for(k = 0; k < numEdgeNode; ++k) 
                    elementEdgeNodeSorted.push_back(elementEdgeNode[i + k]);

                std::pair<int, int> pairTmp(elementTag[i/(numEdgeNode*numNodes)], -1);

                neighbours.push_back(pairTmp);
            }

        }
    
}