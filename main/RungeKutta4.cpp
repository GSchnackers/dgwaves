// Implementation of Runge Kutta of order 4

// Input: - vector u at time step t
//        - time step
//        - derivative of u

// Output: -  vector u at time step t+1

#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"

void RungeKutta4(std::vector<double> & u, const double timeStep, const std::vector<double> & dudt){

    for(std::size_t i=0; i < u.size(); i++){
        k1 = du/dt[i];
        k2 = ;
        k3 = ;
        k4 = ;
        u[i] += timeStep/6 * (k1 + 2*k2 + 2*k3 + k4);
    }
    
}