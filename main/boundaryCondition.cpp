// The functions in this file deal with the boundary conditions.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
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

    std::cout << physicalGroupsTags.size() << std::endl;

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

                    if(physicalName.find("Wall") != std::string::npos)
                        frontierElement.neighbours[j].second = -2;

                    else if(physicalName.find("Constant") != std::string::npos)
                        frontierElement.neighbours[j].second = -3;

                    else if(physicalName.find("Sinusoidal") != std::string::npos)
                        frontierElement.neighbours[j].second = -4;

                    else if(physicalName.find("Out") != std::string::npos)
                        frontierElement.neighbours[j].second = -5;

                    else
                    {
                        gmsh::logger::write("The initial condition is not known.", "error");
                        exit(-1);
                    }
                    
                    
                }
            }
        }

    }

}