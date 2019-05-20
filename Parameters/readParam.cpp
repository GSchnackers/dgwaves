#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <fstream>
#include "structures.hpp"

void readParam(std::string fileName, Simulation & simulation){

    std::ifstream file;
    std::string command;

    simulation.c.resize(3,0);

    if(fileName.find(".wave") == std::string::npos)
    {
        gmsh::logger::write("The parameter file for the simulation is not a .wave. Only .wave are supported.", "error");
        exit(-1);
    }

    file.open(fileName);
    if(!file.is_open())
    {
        gmsh::logger::write("There has been a problem on the opening of the simulation parameter file.", "error");
        exit(-1);
    }

    while(std::getline(file, command, ' '))
    {

        bool cond = true;

        if(!(command.compare("SIMTIME"))) file >> simulation.simTime;
        
        else if(!(command.compare("SIMSTEP"))) file >> simulation.simStep;
        
        else if(!(command.compare("REGISTRATION"))) file >> simulation.registration;
        
        else if(!(command.compare("SOLVER"))) file >> simulation.solver;
        
        else if(!(command.compare("GAUSS"))) 
        {
            std::getline(file, simulation.gaussType);
            cond = false;
        }

        else if(!(command.compare("DEBUG"))) file >> simulation.debug;
    
        else if(!(command.compare("ALPHA"))) file >> simulation.alpha;
        
        else if(!(command.compare("BC")))
        {
            std::getline(file, simulation.boundFileName);
            cond = false;
        }

        else if(!(command.compare("PROP")))
        {
            std::getline(file, simulation.propFileName);
            cond = false;
        }

        else if(!(command.compare("UNUM"))) file >> simulation.uNum;
        
        else if(!(command.compare("C")))
            for(std::size_t i = 0; i < 3; i++)
            {
                file >> simulation.c[i];
                file.get();
            }
        else if(!(command.compare("ERROR"))) file >> simulation.error;

        else if(!(command.compare("NUMTHREADS"))) file >> simulation.numThreads;

        else
            gmsh::logger::write("Unrecognized parameter. Ignored.", "warning");

        if(cond)
            while(1)
            {
                char c = file.get();
                if(c == '\n' || c == EOF) break;
            }

        
    }

    file.close(); 

}