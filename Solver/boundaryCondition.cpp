#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include "structures.hpp"

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
        int physicalNumNodes;

        gmsh::model::mesh::getElementsByType(elemType, physicalElementTag, physicalNodeTags, \
                                             physicalEntityTags[i]);

        physicalNumNodes = physicalNodeTags.size()/physicalElementTag.size();

        for(j = 0; j < frontierElement.elementTag.size(); ++j)
        {
            int check = 0;
            for(k = 0; k < physicalElementTag.size(); ++k)
            {
                int count = 0;
                for(l = 0; l < frontierElement.numNodes; ++l)
                    for(m = 0; m < physicalNumNodes; ++m)
                        if(physicalNodeTags[k * physicalNumNodes + m] == frontierElement.nodeTags[j * frontierElement.numNodes + l])
                            ++count;

                if(count == frontierElement.numNodes)
                {
                    check = 1;
                    break;
                }

                //std::cout << std::endl << count << std::endl;
            }
            

            //std::cout << std::endl;
                
            if(check)
            {
                //std::cout << "Hello" << std::endl;
                //std::cout << frontierElement.elementTag[j] << std::endl;
                while(!std::getline(boundFile, boundCommand, ' ').eof())
                {
                    if(physicalName.find(boundCommand) != std::string::npos)
                    {
                        int mbeg = 0;

                        std::string boundName;

                        std::getline(boundFile, boundName, ' ');

                        if(boundName.find("Sinusoidal") != std::string::npos)
                        {
                            //std::cout << "Fake" << std::endl;
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
                                    int nodeIndex = j * frontierElement.numNodes + l;
                                    int uIndex = frontierElement.neighbours[j].first * \
                                                mainElement.numNodes * 6 + \
                                                frontierElement.nodeCorrespondance[nodeIndex].first * 6 + m;

                                    //std::cout << mainElement.nodeTags[uIndex/6] << " " << frontierElement.nodeTags[nodeIndex] << std::endl;
                                    //std::cout <<  frontierElement.neighbours[j].first << " " << frontierElement.neighbours[j].second << std::endl;
                                    bcParam[uIndex].param1 = par[3 * (m % 3)];
                                    bcParam[uIndex].param2 = par[3 * (m % 3) + 1];
                                    bcParam[uIndex].param3 = par[3 * (m % 3) + 2];
                                    u.boundSign[uIndex] = -2;
                                }

                            frontierElement.neighbours[j].second = -2;
                            std::cout << std::endl;
                            
                        }

                        else if(boundName.find("PerfectCond") != std::string::npos)
                        {

                            for(l = 0; l < frontierElement.numNodes; ++l)
                                for(m = 0; m < 6; ++m)
                                {
                                    int nodeIndex = j * frontierElement.numNodes + l;
                                    int uIndex    = frontierElement.neighbours[j].first * \
                                                    mainElement.numNodes * 6 +\
                                                    frontierElement.nodeCorrespondance[nodeIndex].first * 6 +\
                                                    m;

                                    u.boundSign[uIndex] = -3;
                                }
                            
                            boundFile.get();
                            frontierElement.neighbours[j].second = -3;

                        }

                        else if(boundName.find("Opening") != std::string::npos)
                        {
                            for(l = 0; l < frontierElement.numNodes; ++l)
                                for(m = 0; m < 6; ++m)
                                {
                                    int nodeIndex = j * frontierElement.numNodes + l;
                                    int uIndex    = frontierElement.neighbours[j].first * \
                                                    mainElement.numNodes * 6 +\
                                                    frontierElement.nodeCorrespondance[nodeIndex].first * 6 +\
                                                    m;

                                    u.boundSign[uIndex] = -1;
                                }

                            boundFile.get();
                            frontierElement.neighbours[j].second = -1;
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
                                int nodeIndex = j * frontierElement.numNodes + l;
                                int uIndex    = frontierElement.neighbours[j].first * mainElement.numNodes + \
                                                frontierElement.nodeCorrespondance[nodeIndex].first;

                                bcParam[uIndex].param1 = par[0];
                                bcParam[uIndex].param2 = par[1];
                                bcParam[uIndex].param3 = par[2];

                                u.boundSign[uIndex] = -2;

                            }

                            frontierElement.neighbours[j].second = -2;

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
        if(u.boundSign[i] == -2)
             u.bound[i] = bcParam[i].param1 * sin(bcParam[i].param2 * M_PI * t + bcParam[i].param3);


    for(i = 0; i < u.bound.size() && !t; ++i)
        if(u.boundSign[i] == -3)  u.bound[i] = 0;


}