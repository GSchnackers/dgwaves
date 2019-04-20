/*

This cpp file holds the functions that deals with the boundary conditions imposed on the domain.

-setBoundaryConditions:

    Input: mainElement, a pointer towards the structures of the elements constituying the mesh.
           u, a pointer towards the quantity on which the BC's are applied

    Output: none.

    Working principle: Each point of each element in the domain is associated to a scalar quantity 'u' stored in
                       'u.node'. We associate to each of these scalar quantities a boundary type stored in
                       'u.boundSign'. u.boundSign[i] corresponds to u.node[i] which correspond to 
                       mainElement.nodeTag[i]. 
                       First, the nodes and the name of each physical group are collected. Then, a comparison
                       between the node tag corresponding to u.boundSign[i] and the node tags of the physical 
                       group is performed. If the node tags correspond, then a negative value associated to the 
                       type of the boundary condition applied on the boundary is given to u.boundSign.

- computeBoundaryCondition:

    Input: mainElement and u (same signification as in the function decribed above).
           t, the time of the simulation.

    Output: none

    Working principle: Associates to each node of the boundary its imposed value. The vector containing
                       all information is 'u.bound'.


*/

#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"
#include "structures.h"

void setBoundaryConditions(Element & mainElement, Quantity & u){

    std::size_t i, j, k;
    gmsh::vectorpair physicalGroupsTags;

    // Allow to obtain the dimensions and the tags of the physical groups.
    // It is assumed here that all physical groups have the same dimension (1 or 2), representing
    // where the BC's have to be applied.

    gmsh::model::getPhysicalGroups(physicalGroupsTags, mainElement.dim - 1);

    // Run through all physical groups of dimension 1 in 2D and 2 in 3D.
    for(i = 0; i < physicalGroupsTags.size(); ++i)
    {

        std::vector<int> physicalNodeTags;
        std::vector<double> coords;

        std::string physicalName;

        gmsh::model::mesh::getNodesForPhysicalGroup(physicalGroupsTags[i].first, physicalGroupsTags[i].second,\
                                                    physicalNodeTags, coords);

        gmsh::model::getPhysicalName(physicalGroupsTags[i].first, physicalGroupsTags[i].second, physicalName);

        // Run through each edge.
        for(j = 0; j < physicalNodeTags.size(); ++j)
            for(k = 0; k < mainElement.nodeTags.size(); ++k)
                if(mainElement.nodeTags[k] == physicalNodeTags[j])
                {
                    if(physicalName.find("Sinusoidal") != std::string::npos)
                        u.boundSign[k] = -2;

                    else if(physicalName.find("Constant") != std::string::npos)
                        u.boundSign[k] = -3;

                }

    }

}


void computeBoundaryCondition(Quantity & u, const double t){

    std::size_t i;

    for(i = 0; i < u.bound.size(); ++i)
    {

        if(u.boundSign[i] == -2)
             u.bound[i] = sin(2 * M_PI * t/0.5);
        
    
        else if(u.boundSign[i] == -3)
             u.bound[i] = 1;

    }


}