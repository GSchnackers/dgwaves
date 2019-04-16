#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"
#include "structures.h"

// This function set the specific type of boundary condition applied to the specific place of the frontier.
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

                    //std::cout << mainElement.nodeTags[k] << std::endl;

                }

    }

}


void computeBoundaryCondition(const Element & mainElement, Quantity & u, const double t){

    std::size_t i;

    for(i = 0; i < mainElement.nodeTags.size(); ++i)
    {

        if(u.boundSign[i] == -2)
             u.bound[i] = sin(2 * M_PI * t/0.5);
        
    
        else if(u.boundSign[i] == -3)
             u.bound[i] = 1;

    }


}