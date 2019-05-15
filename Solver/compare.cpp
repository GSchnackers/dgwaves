#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "solver.hpp"
#include "structures.hpp"

void compare(std::vector<double> & error, std::vector<double> & errorNodes, Quantity & u,\
             const std::vector<double> & coordinates, const Element & mainElement,\
             const Simulation & simulation, const double mytime, const Properties & matProp){
    
    int numNodes = u.node.size()/simulation.uNum;
    std::vector<double> uAnal(numNodes*simulation.uNum);
    double w, wc, beta, c, muR, epsR;
    double a = 1;
    int m = 1;
    double f = 2;
    
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
            muR = matProp.relPermeability.node[i];
            epsR = matProp.relPermittivity.node[i];
            c = 1/sqrt(epsR*muR);
            wc = c*m*M_PI/a;
            w = 2*M_PI*f;
            beta = sqrt(w*w-wc*wc)/c;
            if(coordinates[3*i] < mytime){
                // Dimensionless amplitude is taken = 1 (to change in the code if any change in the BC)
                uAnal[simulation.uNum*i] = 0; // Ex
                uAnal[simulation.uNum*i + 1] = 0; // Ey
                uAnal[simulation.uNum*i + 2] = sin(w*mytime-beta*coordinates[3*i])*sin(m*M_PI*coordinates[3*i+1]/a); // Ez
                uAnal[simulation.uNum*i + 3] = -(m*M_PI/(a*w*muR))*cos(w*mytime-beta*coordinates[3*i])*cos(m*M_PI*coordinates[3*i+1]/a); // Hx
                uAnal[simulation.uNum*i + 4] = -(beta/(w*muR))*sin(w*mytime-beta*coordinates[3*i])*sin(m*M_PI*coordinates[3*i+1]/a); // Hy
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

        for(std::size_t j = 0; j < simulation.uNum ; j++){
            errorNodes[simulation.uNum*i + j] = uAnal[simulation.uNum*i + j] - u.node[simulation.uNum*i + j];
            error[j] += errorNodes[simulation.uNum*i + j]*errorNodes[simulation.uNum*i + j]/numNodes;
        }
    }
}