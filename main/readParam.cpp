#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <fstream>
#include "functions.h"
#include "structures.h"

void readParam(std::string fileName, double & simTime, double & incrementation, int & registration,
               int & solvType, std::string & gaussType, int & meshDim, int & debug, double & alpha){

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
             file >> simTime;
             file.get();
        }
        
        else if(!(command.compare("SIMSTEP")))
        {
            file >> incrementation;
            file.get();
        }

        else if(!(command.compare("REGISTRATION")))
        {
            file >> registration;
            file.get();
        }

        else if(!(command.compare("SOLVER")))
        {
            file >> solvType;
            file.get();
        }

        else if(!(command.compare("GAUSS"))) std::getline(file, gaussType);

        else if(!(command.compare("MESHDIM")))
        { 
            file >> meshDim;
            file.get();
        }

        else if(!(command.compare("DEBUG")))
        { 
            file >> debug;
            file.get();
        }

        else if(!(command.compare("ALPHA")))
        { 
            file >> alpha;
            file.get();
        }

        else
            gmsh::logger::write("Unrecognized parameter. Ignored.", "warning");
        
    }

    file.close(); 

}