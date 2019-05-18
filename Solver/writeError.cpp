#include <cstdio>
#include <iostream>
#include <fstream>
#include <gmsh.h>
#include "solver.hpp"
#include "structures.hpp"

void writeError(const std::vector<double> & error, const Simulation & simulation){

    std::ofstream fichier("error.txt", std::ios::out | std::ios::trunc);
    
    if(fichier) {
        if(simulation.uNum == 1)
            fichier << "time error Power" << std::endl;
        else if(simulation.uNum == 6)
            fichier << "time errorEx errorEy errorEz errorHx errorHy errorHz Power" << std::endl;
        
        for(std::size_t i = 0; i < error.size()/(simulation.uNum+1); i++){
            fichier << std::to_string(i*simulation.simStep);
            for(std::size_t j = 0; j < simulation.uNum+1; j++)
                fichier << " " << std::to_string(error[i*(simulation.uNum+1) + j]);
            fichier << std::endl;
        }

        fichier.close();
    }
    else
        std::cerr << "Error to opening!" << std::endl;
}