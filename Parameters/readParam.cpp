#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <fstream>
#include "structures.hpp"

void readParam(std::string fileName, Simulation & simulation){

    std::ifstream file;
    std::string command;

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

        if(!(command.compare("SIMTIME")))
        {
             file >> simulation.simTime;
             file.get();
        }
        
        else if(!(command.compare("SIMSTEP")))
        {
            file >> simulation.simStep;
            file.get();
        }

        else if(!(command.compare("REGISTRATION")))
        {
            file >> simulation.registration;
            file.get();
        }

        else if(!(command.compare("SOLVER")))
        {
            file >> simulation.solver;
            file.get();
        }

        else if(!(command.compare("GAUSS"))) std::getline(file, simulation.gaussType);

        else if(!(command.compare("DEBUG")))
        { 
            file >> simulation.debug;
            file.get();
        }

        else if(!(command.compare("ALPHA")))
        { 
            file >> simulation.alpha;
            file.get();
        }

        else if(!(command.compare("E0")))
        { 
            file >> simulation.E0;
            file.get();
        }

        else if(!(command.compare("L")))
        { 
            file >> simulation.L;
            file.get();
        }

        else if(!(command.compare("BC"))) std::getline(file, simulation.boundFileName);

        else if(!(command.compare("PROP"))) std::getline(file, simulation.propFileName);

        else if(!(command.compare("UNUM")))
        { 
            file >> simulation.uNum;
            file.get();
        }

        else
            gmsh::logger::write("Unrecognized parameter. Ignored.", "warning");
        
    }

    file.close(); 

}