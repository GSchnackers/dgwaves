#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

void compare(double & error, std::vector<double> & errorNodes, const std::vector<double> coefF,\
            const std::vector<double> & coordinates, std::vector<double> & u, const double mytime){
                
    int numNodes = u.size();
    std::vector<double> uAnal(numNodes);
    double T = 0.5;
    double w = 2*M_PI/T;
    error = 0;

    for(std::size_t i = 0; i < numNodes; i++){
        if(coordinates[3*i] < coefF[0]*mytime)
            uAnal[i] = sin(w*mytime- coordinates[3*i]*w/coefF[0]);
        else
            uAnal[i] = 0;
        
        errorNodes[i] = uAnal[i] - u[i];
        error += errorNodes[i]*errorNodes[i]/numNodes;
    }
}