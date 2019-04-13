// The functions in this file deal with the boundary conditions.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"
#include "structures.h"

// This function set the specific type of boundary condition applied to the specific place of the frontier.
void setBoundaryConditions(Element & mainElement, Quantity & u){

    std::size_t i, j, k, l;
    gmsh::vectorpair physicalGroupsTags;

    // Allow to obtain the dimensions and the tags of the physical groups.
    // It is assumed here that all physical groups have the same dimension (1 or 2), representing
    // where the BC's have to be applied.

    gmsh::model::getPhysicalGroups(physicalGroupsTags, mainElement.dim - 1);

    for(i = 0; i < physicalGroupsTags.size(); ++i)
    {

        std::vector<int> physicalNodeTags;
        std::vector<double> coords;

        std::string physicalName;

        gmsh::model::mesh::getNodesForPhysicalGroup(physicalGroupsTags[i].first, physicalGroupsTags[i].second,\
                                                    physicalNodeTags, coords);

        gmsh::model::getPhysicalName(physicalGroupsTags[i].first, physicalGroupsTags[i].second, physicalName);

        for(j = 0; j < mainElement.nodeTags.size(); ++j)
        {
            for(k = 0; k < physicalNodeTags.size(); ++k)
                if(physicalNodeTags[k] == mainElement.nodeTags[j])
                {
                
                    if(physicalName.find("Sinusoidal") != std::string::npos)
                        u.bound[j] = -1;
                    if(physicalName.find("Constant") != std::string::npos)
                        u.bound[j] = -2;
                    
                }
            
        }

    }

}

// At first, very simple boundary conditions. They are applied at the gauss points of the boundaries of the
// domain, since the fluxes are applied there. The effects of those boundary conditions are only felt at the
// Gauss points of the frontier element.
void computeBoundaryCondition(const Element & mainElement, Quantity & u, const double t){

    std::size_t i, j;

    for(i = 0; i < mainElement.nodeTags.size(); ++i)
    {
        if(u.bound[i] == -1)
             u.node[i] = sin(t);
        if(u.bound[i] == -2)
             u.node[i] = 1;        
    }


}