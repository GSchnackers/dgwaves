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

void boundAssign(Element & frontierElement, const Element & mainElement, const std::vector<int> & physicalEntityTags, \
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
                        if(physicalNodeTags[j * frontierElement.numNodes + l] == \
                           frontierElement.nodeTags[k * frontierElement.numNodes + m]) 
                            ++count;

                
                if(count == frontierElement.numNodes)
                {
                    std::cout << frontierElement.elementTag[k] << std::endl;
                    while(!std::getline(boundFile, boundCommand, ' ').eof())
                    {
                        if(physicalName.find(boundCommand) != std::string::npos)
                        {
                            int mbeg = 0;

                            std::string boundName;

                            std::getline(boundFile, boundName, ' ');

                            if(boundName.find("Sinusoidal") != std::string::npos)
                            {
                                std::vector<double> par(9);

                                for(l = 0; l < par.size(); ++l)
                                {
                                    boundFile >> par[l];
                                    boundFile.get();
                                }

                                if(boundName.find("SinusoidalH") != std::string::npos) mbeg = 3;

                                for(l = 0; l < frontierElement.numNodes; ++l)
                                    for(m = mbeg; m < mbeg + 3; ++m)
                                    {
                                        int nodeIndex = k * frontierElement.numNodes + l;
                                        int uIndex = frontierElement.neighbours[k].first * \
                                                     mainElement.numNodes * 6 +\
                                                     frontierElement.nodeCorrespondance[nodeIndex].first * 6 + m;

                                        bcParam[uIndex].param1 = par[3 * (m % 3)];
                                        bcParam[uIndex].param2 = par[3 * (m % 3) + 1];
                                        bcParam[uIndex].param3 = par[3 * (m % 3) + 2];
                                        u.boundSign[uIndex] = -2;
                                    }

                                frontierElement.neighbours[k].second = -2;
                                
                            }

                            else if(boundName.find("PerfectCond") != std::string::npos)
                            {
                                for(l = 0; l < frontierElement.numNodes; ++l)
                                    for(m = 0; m < 6; ++m)
                                    {
                                        int nodeIndex = k * frontierElement.numNodes + l;
                                        int uIndex = frontierElement.neighbours[k].first * \
                                                     mainElement.numNodes * 6 +\
                                                     frontierElement.nodeCorrespondance[nodeIndex].first * 6 +\
                                                     m;

                                        u.boundSign[uIndex] = -3;
                                    }
                                
                                boundFile.get();
                                frontierElement.neighbours[k].second = -3;
                            }

                            else if(boundName.find("Opening") != std::string::npos)
                            {
                                for(l = 0; l < frontierElement.numNodes; ++l)
                                    for(m = 0; m < 6; ++m)
                                    {
                                        int nodeIndex = k * frontierElement.numNodes + l;
                                        int uIndex    = frontierElement.neighbours[k].first * \
                                                        mainElement.numNodes * 6 +\
                                                        frontierElement.nodeCorrespondance[nodeIndex].first * 6 +\
                                                        m;

                                        u.boundSign[uIndex] = -4;
                                    }

                                boundFile.get();
                                frontierElement.neighbours[k].second = -4;
                            }

                            else if(boundName.find("Sine") != std::string::npos)
                            {
                                std::vector<int> par(3);

                                for(l = 0; l < 3; ++l)
                                {
                                    boundFile >> par[l];
                                    boundFile.get();
                                }

                                for(l = 0; l < frontierElement.numNodes; ++l)
                                {
                                    int nodeIndex = k * frontierElement.numNodes + l;
                                    int uIndex    = frontierElement.neighbours[k].first * mainElement.numNodes + \
                                                    frontierElement.nodeCorrespondance[nodeIndex].first;

                                    bcParam[uIndex].param1 = par[0];
                                    bcParam[uIndex].param2 = par[1];
                                    bcParam[uIndex].param3 = par[2];

                                    u.boundSign[uIndex] = -2;

                                }

                                frontierElement.neighbours[k].second = -2;
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

void setBoundaryCondition(Element & frontierElement, const Element & mainElement,\
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
        if(physicalGroups.dimTags[i].first == frontierElement.dim)
            for(j = 0; j < physicalGroups.elemType[i][0].size(); ++j)
                boundAssign(frontierElement, mainElement, physicalGroups.entityTags[i], boundaryFile, \
                            physicalGroups.name[i], physicalGroups.elemType[i][0][j], frontierElement.dim, \
                            bcParam, u);

    boundaryFile.close();

}


void computeBoundaryCondition(Quantity & u, const double t, const std::vector<Parameter> & bcParam){

    std::size_t i;

    for(i = 0; i < u.bound.size(); ++i)
    {
        
        if(u.boundSign[i] == -2)
             u.bound[i] = bcParam[i].param1 * sin(bcParam[i].param2 * M_PI * t + bcParam[i].param3);
        
    
        // In the case of a perfect conductor, the electric and magnetic field are assumed to be 0.
        else if(u.boundSign[i] == -3)
             u.bound[i] = 0;

    }


}