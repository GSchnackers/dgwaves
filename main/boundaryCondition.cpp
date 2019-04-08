// The functions in this file deal with the boundary conditions.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"
#include "structures.h"

// This function set the specific type of boundary condition applied to the specific place of the frontier.
void setBoundaryConditions(Element & frontierElement){

    std::size_t i, j, k, l;
    gmsh::vectorpair physicalGroupsTags;

    // Allow to obtain the dimensions and the tags of the physical groups.
    // It is assumed here that all physical groups have the same dimension (1 or 2), representing
    // where the BC's have to be applied.

    gmsh::model::getPhysicalGroups(physicalGroupsTags);

    for(i = 0; i < physicalGroupsTags.size(); ++i)
    {

        std::vector<int> physicalNodeTags;
        std::vector<double> coords;

        std::string physicalName;

        gmsh::model::mesh::getNodesForPhysicalGroup(physicalGroupsTags[i].first, physicalGroupsTags[i].second,\
                                                    physicalNodeTags, coords);

        gmsh::model::getPhysicalName(physicalGroupsTags[i].first, physicalGroupsTags[i].second, physicalName);

        for(j = 0; j < frontierElement.elementTag.size(); ++j)
        {

            int count = 0;

            for(k = 0; k < frontierElement.numNodes; ++k)
            {

                for(l = 0; l < physicalNodeTags.size(); ++l)
                {
                    int frontIndex = j * frontierElement.numNodes + k;
                    
                    if(frontierElement.nodeTags[frontIndex] == physicalNodeTags[l]) ++count;
                }

                if(count == frontierElement.numNodes)
                {

                    if(physicalName.find("Output") != std::string::npos)
                        frontierElement.neighbours[j].second = -1;

                    else if(physicalName.find("Wall") != std::string::npos)
                        frontierElement.neighbours[j].second = -2;

                    else if(physicalName.find("Constant") != std::string::npos)
                        frontierElement.neighbours[j].second = -3;

                    else if(physicalName.find("Sinusoidal") != std::string::npos)
                        frontierElement.neighbours[j].second = -4;

                    else
                    {
                        gmsh::logger::write("The initial condition is not known.", "error");
                        std::cout << std::string::npos << std::endl;
                        exit(-1);
                    }
                    
                }
            }
        }

    }

}

// At first, very simple boundary conditions. They are applied at the gauss points of the boundaries of the
// domain, since the fluxes are appleid there.
void computeBoundaryCondition(const Element & frontierElement, Quantity & u, const double t){

    std::size_t i, j;

    u.numGp.resize(frontierElement.elementTag.size() * frontierElement.numGp);

    for(i = 0; i < frontierElement.elementTag.size(); ++i)
        for(j = 0; j < frontierElement.numGp; ++j)
        {
            int index = i * frontierElement.numGp + j;
            if(t == 0) // The fixed boundary condition with time do not need to be revaluated farther.
                switch (frontierElement.neighbours[i/frontierElement.numNodes].second)
                {
                    case -2:
                        u.numGp[index].second = 0;
                        break;

                    case -3:
                        u.numGp[index].second = 1;
                        break;
                }

            if(frontierElement.neighbours[i].second == -4)
                u.numGp[index].second = sin(t);
        }

        

}