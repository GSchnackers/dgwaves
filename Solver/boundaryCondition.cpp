#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include "structures.hpp"

void boundAssign(Element & frontierElement, const std::vector<int> & physicalEntityTags, \
                std::fstream & boundFile, const std::string & physicalName,\
                int elemType, int elemDim, Quantity & u){

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
            // This loop check wether a frontier element is in the physical group or not.
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

            }
            
            //  If the element is in a physical group, a BC type is assigned to it thanks to the .bc file.
            if(check)
            {
                while(!std::getline(boundFile, boundCommand, ' ').eof())
                {
                    if(boundCommand.compare(physicalName) != std::string::npos)
                    {

                        std::string boundName;
                        int l = 0;
                        char tmp = ' ';

                        std::getline(boundFile, boundName, ' ');

                        if(!boundName.compare("SinusoidalE"))
                            frontierElement.neighbours[j].second = SINUS_E;
                        
                        else if(!boundName.compare("SinusoidalH"))
                            frontierElement.neighbours[j].second = SINUS_H;

                        else if(!boundName.compare("PerfectCond"))
                            frontierElement.neighbours[j].second = PERFECTCOND;

                        else if(!boundName.compare("Sine"))
                            frontierElement.neighbours[j].second = SINE;
                        
                        else if(!boundName.compare("TE2D"))
                            frontierElement.neighbours[j].second = TE2D;

                        else if(!boundName.compare("TE3D"))
                            frontierElement.neighbours[j].second = TE3D;

                        boundFile >> frontierElement.bcParam[j * 9 + l];
                        ++l;

                        while(1){
                        
                            char c;
                            c = boundFile.get();
                            if(c == EOF || c == '\n') break;
                            boundFile >> frontierElement.bcParam[j * 9 + l];
                            ++l;
                            
                        }

                        break;

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

void setBoundaryCondition(Element & frontierElement, const Simulation & simulation, \
                          const PhysicalGroups & physicalGroups, Quantity & u){

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
                boundAssign(frontierElement, physicalGroups.entityTags[i], boundaryFile, \
                            physicalGroups.name[i], physicalGroups.elemType[i][0][j], frontierElement.dim, u);

    boundaryFile.close();

}