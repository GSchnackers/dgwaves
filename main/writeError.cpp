#include <cstdio>
#include <iostream>
#include <fstream>
#include <gmsh.h>
#include "functions.h"

void writeError(std::vector<double> & error, const double timeStep) {

    std::ofstream fichier("error.txt", std::ios::out | std::ios::trunc);
    
    if(fichier) {
        fichier << "time error" << std::endl;
        for(std::size_t i = 0; i < error.size(); i++)
            fichier << std::to_string(i*timeStep) << " " << std::to_string(error[i]) << std::endl;

        fichier.close();
    }
    else
        std::cerr << "Error to opening!" << std::endl;
}