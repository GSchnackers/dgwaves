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

void boundAssign(const Element & frontierElement, const Element & mainElement, const std::vector<int> & physicalEntityTags, \
                std::fstream & boundFile, const std::string & physicalName,\
                int elemType, int elemDim, std::vector<Parameter> & bcParam, Quantity & u){

    std::size_t i, j, k, l, m;

    std::vector<int> binInt;
    std::string bin;
    std::string boundCommand;

    // For each physical entity tags, 
    for(i = 0; i < physicalEntityTags.size(); ++i)
    {
        std::vector<int> physicalElementTag;
        std::vector<int> physicalNodeTags;
        std::vector<double> bin1, bin2;

        gmsh::model::mesh::getElementsByType(elemType, physicalElementTag, physicalNodeTags, \
                                             physicalEntityTags[i]);
        
        for(j = 0; j < physicalElementTag.size(); ++j)
            for(k = 0; k < frontierElement.elementTag.size(); ++k)
            {
                int count = 0;
                for(l = 0; l < frontierElement.numNodes; ++l)
                    for(m = 0; m < frontierElement.numNodes; ++m)
                        if(physicalNodeTags[j * frontierElement.numNodes + m] == \
                           frontierElement.nodeTags[k * frontierElement.numNodes + l]) 
                            ++count;

                std::cout << count << std::endl;
                
                if(count == frontierElement.numNodes)
                {
                    //std::cout << "Hello" << std::endl;
                    while(!std::getline(boundFile, boundCommand, ' ').eof())
                    {
                        if(physicalName.find(boundCommand) != std::string::npos)
                        {
                            int mbeg = 0, type = -1;

                            std::string boundName;

                            std::getline(boundFile, boundName, ' ');

                            std::cout << boundName << std::endl;

                            if(boundName.find("Sinusoidal") != std::string::npos)
                            {
                                std::vector<double> par(9);

                                for(l = 0; l < par.size(); ++l)
                                {
                                    boundFile >> par[l];
                                    boundFile.get();
                                }

                                if(boundName.find("SinusoidalE") != std::string::npos) mbeg = 0;

                                else if(boundName.find("SinusoidalH") != std::string::npos) mbeg = 3;

                                for(l = 0; l < frontierElement.numNodes; ++l)
                                    for(m = mbeg; m < mbeg + 3; ++m)
                                    {
                                        int uIndex = frontierElement.neighbours[k].first * mainElement.numNodes * 6 +\
                                                     frontierElement.nodeCorrespondance[k * \
                                                                 frontierElement.numNodes + l].first * 6 +\
                                                     m;

                                        bcParam[uIndex].param1 = par[3 * (m % 3)];
                                        bcParam[uIndex].param2 = par[3 * (m % 3) + 1];
                                        bcParam[uIndex].param3 = par[3 * (m % 3) + 2];
                                        u.boundSign[uIndex] = -2;
                                    }
                                
                            }

                        }

                        else
                            std::getline(boundFile, bin);
                    }

                    boundFile.clear();
                    boundFile.seekg(0, std::ios::beg);

                }
            }                
        
    }

}

void setBoundaryCondition(const Element & frontierElement, const Element & mainElement,\
                          const Simulation & simulation, const PhysicalGroups & physicalGroups, \
                          Quantity & u, std::vector<Parameter> & bcParam){

    std::size_t i, j;

    std::fstream boundaryFile;

    if(simulation.boundFileName.find(".bc") == std::string::npos)
    {
        gmsh::logger::write("Unsupported boundary condition file extension.", "error");
        exit(-1);
    }

    boundaryFile.open(simulation.boundFileName);

    if(!boundaryFile.is_open())
    {
        gmsh::logger::write("Could not open the boundary condition file.", "error");
        exit(-1);
    }

    
    for(i = 0; i < physicalGroups.dimTags.size(); ++i)
        if(physicalGroups.dimTags[i].second == frontierElement.dim)
            for(j = 0; j < physicalGroups.elemType[i][0].size(); ++j)
            {  
                boundAssign(frontierElement, mainElement, physicalGroups.entityTags[i], boundaryFile, \
                            physicalGroups.name[i], physicalGroups.elemType[i][0][j], frontierElement.dim, \
                            bcParam, u);
                //std::cout << physicalGroups.elemType[i][0][j] << std::endl;
            }

    boundaryFile.close();

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