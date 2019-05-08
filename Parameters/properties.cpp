#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include <fstream>
#include "structures.hpp"

void propAssign(const Element & element, const std::vector<int> & physicalEntityTags, \
                std::fstream & propFile, const std::string & physicalName, const int elemDim,\
                Properties & matProp){

    std::size_t i, j, k, l;

    std::vector<int> binInt;
    std::string bin;
    std::string propCommand;

    // For each physical entity tags, 
    for(i = 0; i < physicalEntityTags.size(); ++i)
    {
        std::vector<int> physicalElementTag;

        gmsh::model::mesh::getElementsByType(elemDim, physicalElementTag, binInt, \
                                             physicalEntityTags[i]);

        for(j = 0; j < physicalElementTag.size(); ++j)
            for(k = 0; k < element.elementTag.size(); ++k)
                if(physicalElementTag[j] == element.elementTag[k])
                {
                    while(!std::getline(propFile, propCommand, ' ').eof())
                    {
                        if(physicalName.find(propCommand) != std::string::npos)
                        {
                            double relPermitt, relPermea, conduct;
                            propFile >> relPermitt;
                            propFile.get();
                            propFile >> relPermea;
                            propFile.get();
                            propFile >> conduct;
                            propFile.get();

                            for(l = 0; l < element.numNodes; ++l)
                            {

                                matProp.relPermittivity.bound[k * element.numNodes + l] = \
                                matProp.relPermittivity.node[k * element.numNodes + l] = relPermitt;

                                matProp.relPermeability.bound[k * element.numNodes + l] = \
                                matProp.relPermeability.node[k * element.numNodes + l] = relPermea;

                                matProp.conductivity.bound[k * element.numNodes + l] = \
                                matProp.conductivity.node[k * element.numNodes + l] = conduct;

                            }


                        }

                        else
                            std::getline(propFile, bin);


                    }

                    propFile.clear();
                    propFile.seekg(0, std::ios::beg);

                }
                            
        
    }

}

void setProperties(const Element & mainElement, const Element & frontierElement, const Simulation & simulation, \
                   const PhysicalGroups & physicalGroups, Properties & matProp){


    std::size_t i, j;
    
    std::fstream propFile;

    propFile.open(simulation.propFileName);
    if(simulation.propFileName.find(".prop") == std::string::npos)
    {
        gmsh::logger::write("Unsupported properties file extension.", "error");
        exit(-1);
    }

    if(!propFile.is_open())
    {
        gmsh::logger::write("Could not open the property file.", "error");
        exit(-1);
    }

    for(i = 0; i < physicalGroups.dimTags.size(); ++i)
    {

        // Assigns the properties of the main elements.
        for(j = 0; j < physicalGroups.elemType[i][0].size(); ++j)
        {
            propAssign(mainElement, physicalGroups.entityTags[i], propFile, physicalGroups.name[i],\
                    physicalGroups.elemType[i][0][j], matProp);
        
        // Assigns the properties of the frontier elements.
            propAssign(frontierElement, physicalGroups.entityTags[i], propFile,  physicalGroups.name[i], \
                    physicalGroups.elemType[i][0][j], matProp);
        }

    }

    propFile.close();
}