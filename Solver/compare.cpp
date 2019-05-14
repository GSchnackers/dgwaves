#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "solver.hpp"
#include "structures.hpp"

void compare(double & error, std::vector<double> & errorNodes, Quantity & u,\
             const std::vector<double> & coordinates, const Element & mainElement,\
             const Simulation & simulation, const double mytime){
    
    int numNodes = u.node.size()/simulation.uNum;
    std::vector<double> uAnal(numNodes*simulation.uNum);
    double w;
    error = 0;
    
    for(std::size_t i = 0; i < numNodes; i++){
        if(simulation.uNum == 1){
            if(coordinates[3*i] < simulation.c[0]*mytime){
                //w = bcParam[0].param2 * M_PI;
                //uAnal[i] = bcParam[0].param1 * \
                sin(w * (mytime - coordinates[3*i]/simulation.c[0]) + bcParam[0].param3);
                uAnal[i] = 0;
            }
            else
                uAnal[i] = 0;
        } else if(simulation.uNum == 6) {
            if(0){
                uAnal[simulation.uNum*i] = 0; // Ex
                uAnal[simulation.uNum*i + 1] = 0; // Ey
                uAnal[simulation.uNum*i + 2] = 0; // Ez
                uAnal[simulation.uNum*i + 3] = 0; // Hx
                uAnal[simulation.uNum*i + 4] = 0; // Hy
                uAnal[simulation.uNum*i + 5] = 0; // Hz
            }
            else{
                uAnal[simulation.uNum*i] = 0; // Ex
                uAnal[simulation.uNum*i + 1] = 0; // Ey
                uAnal[simulation.uNum*i + 2] = 0; // Ez
                uAnal[simulation.uNum*i + 3] = 0; // Hx
                uAnal[simulation.uNum*i + 4] = 0; // Hy
                uAnal[simulation.uNum*i + 5] = 0; // Hz
            }
        }

        errorNodes[i] = uAnal[i] - u.node[i];
        error += errorNodes[i]*errorNodes[i]/numNodes;
    }
}