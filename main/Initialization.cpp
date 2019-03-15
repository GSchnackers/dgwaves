#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void Initialization(std::vector<struct Entity> & geometry){

    std::size_t i, j; // Index variables
    gmsh::vectorpair entities2D; // Vector pair used to contain the dimension and tags of entities.
    gmsh::model::getEntities(entities2D, 2); // Load the 2D entities tag in the memory.
    std::cout << entities2D.size();
    
    std::vector<struct Entity> tmpGeo(entities2D.size()); // Temporary vector containing all 2D entities.

    // Initialization of the 2D elements.

    for(i = 0; i < tmpGeo.size(); ++i){

        tmpGeo[i].entityTag = entities2D[i].second;

        gmsh::model::mesh::getElementTypes(tmpGeo[i].entityElementProperties2D.elementType, 2);
        if(tmpGeo[i].entityElementProperties2D.elementType.size() > 1){

            gmsh::logger::write("No hybrid implementation is authorized.", "error");
            exit(-1);

        }

        // Gets 
        gmsh::model::mesh::getElementsByType(tmpGeo[i].entityElementProperties2D.elementType[0],\
                                            tmpGeo[i].elementTag2D, tmpGeo[i].nodeTags,\
                                            tmpGeo[i].entityTag);
                                            
        gmsh::model::mesh::getElementEdgeNodes(tmpGeo[i].entityElementProperties2D.elementType[0],\
                                               tmpGeo[i].elementEdgeNode2D, tmpGeo[i].entityTag);

        //edges(element2D.numNodes, tmpGeo[i]);

    }

    geometry = tmpGeo;

     // Gets the type of all elements of dimension 2 of all tags

}

// REM: can be easily modified to comprehend hybrid elements.