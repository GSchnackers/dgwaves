// Implementation of the forward euler method

// Input: - vector u at time step t
//        - timestep
//        - derivative of u

// Output: -  vector u at time step t+1

#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"

void Forward_Euler_method(std::vector<double> & u, double timestep, std::vector<double> dudt){

    for(std::size_t i=0; i < u.size();i++)
        u[i]=u[i]+dudt[i]*timestep;
    
}