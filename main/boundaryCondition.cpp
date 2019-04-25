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
#include <fstream>
#include "functions.h"
#include "structures.h"

void setBoundaryConditions(const Element & mainElement, const std::string & boundaryFileName,\
                           const std::string & propFileName, Quantity & relPermittivity,\
                           Quantity & relPermeability, Quantity & conductivity, Quantity & u,\
                           std::vector<Parameter> & bcParam){

    std::size_t i, j, k;
    gmsh::vectorpair physicalGroupsTags;

    std::fstream boundaryFile;
    std::fstream propFile;

    std::string boundCommand;
    std::string propCommand;

    if(boundaryFileName.find(".bc") == std::string::npos || propFileName.find(".prop") == std::string::npos)
    {
        gmsh::logger::write("Unsupported file extension.", "error");
        exit(-1);
    }

    boundaryFile.open(boundaryFileName);

    if(!boundaryFile.is_open())
    {
        gmsh::logger::write("Could not open the boundary condition file.", "error");
        exit(-1);
    }

    propFile.open(propFileName);
    if(!propFile.is_open())
    {
        gmsh::logger::write("Could not open the property file.", "error");
        boundaryFile.close();
        exit(-1);
    }

    gmsh::model::getPhysicalGroups(physicalGroupsTags);

    while(std::getline(boundaryFile, boundCommand, ' '))
        while(std::getline(propFile, propCommand, ' ')) // Run through all physical groups of dimension 1 in 2D and 2 in 3D.
            for(i = 0; i < physicalGroupsTags.size(); ++i)
            {

                std::vector<int> physicalNodeTags;
                std::vector<double> coords;

                std::string physicalName;

                gmsh::model::mesh::getNodesForPhysicalGroup(physicalGroupsTags[i].first,\
                                                            physicalGroupsTags[i].second,\
                                                            physicalNodeTags, coords);

                gmsh::model::getPhysicalName(physicalGroupsTags[i].first, physicalGroupsTags[i].second,\
                                             physicalName);

                // Run through each edge.
                for(j = 0; j < physicalNodeTags.size(); ++j)
                    for(k = 0; k < mainElement.nodeTags.size(); ++k)
                        if(mainElement.nodeTags[k] == physicalNodeTags[j])
                        {
                            if(!physicalName.compare(boundCommand))
                            {
                                std::string boundName;

                                std::getline(boundaryFile, boundName, ' ');

                                if(!boundName.compare("Sinusoidal"))
                                {
                                    u.boundSign[k] = -2;
                                    boundaryFile >> bcParam[k].param1;
                                    boundaryFile.get();
                                    boundaryFile >> bcParam[k].param2;
                                    boundaryFile.get();
                                    boundaryFile >> bcParam[k].param3;
                                    
                                }

                                else if(!boundName.find("Constant"))
                                {
                                    u.boundSign[k] = -3;
                                    boundaryFile >> bcParam[k].param1;
                                }

                                else
                                    gmsh::logger::write("Unknown BC type. Ignored.", "warning");
                                
                                
                            }

                            else if(!physicalName.compare(propCommand))
                            {
                                propFile >> relPermittivity.node[k];
                                propFile.get();
                                propFile >> relPermeability.node[k];
                                propFile.get();
                                propFile >> conductivity.node[k];
                            }

                            boundaryFile.get();

                        }

            }
        
    
    boundaryFile.close();
    propFile.close();

}


void computeBoundaryCondition(Quantity & u, const double t, const std::vector<Parameter> & bcParam){

    std::size_t i;

    for(i = 0; i < u.bound.size(); ++i)
    {

        if(u.boundSign[i] == -2)
             u.bound[i] = bcParam[i].param1 * sin(bcParam[i].param2 * M_PI * t + bcParam[i].param3);
        
    
        else if(u.boundSign[i] == -3)
             u.bound[i] = bcParam[i].param1;

    }


}