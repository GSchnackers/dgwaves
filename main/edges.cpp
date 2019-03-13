#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"
#include "structures.h"

void edges(const int numNodes, struct Entity & entity)
{
    std::size_t i, j;

    entity.edges.push_back(entity.elementEdgeNodes[0]);
    entity.edges.push_back(entity.elementEdgeNodes[1]);

    for (i = 0; i < entity.elementEdgeNodes.size(); i += 2)

        for (j = 0; j < entity.edges.size(); j += 2)
        {
            // Check if edges is already in entity.edges in the same or opposite direction.
            bool condition = (entity.elementEdgeNodes[i] == entity.edges[j] \
                             && entity.elementEdgeNodes[i + 1] == entity.edges[j + 1]) ||\
                            (entity.elementEdgeNodes[i] == entity.edges[j + 1] \
                            && entity.elementEdgeNodes[i + 1] == entity.edges[j]);

            if(condition){

                entity.neighbours[j + 1] = entity.elementTags[j/(2*numNodes)];  // ok seulement pour T3 car si T6 : le noeud du milieu du coté n'apparait qu'une seule fois
                break;

            }

            // If the edge is not already in entity.edges, we add it to the sorting node.
            if(j + 2 == entity.edges.size())
            {
                entity.edges.push_back(entity.elementEdgeNodes[i]);
                entity.edges.push_back(entity.elementEdgeNodes[i+1]);
                entity.neighbours.push_back(entity.elementTags[i/(2*numNodes)]);
                entity.neighbours.push_back(-1);
            }

        }
    
}