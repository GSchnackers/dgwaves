#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void boundaryConditions(std::vector<double> & coord,double mytime, double & value){
    
    //value = 1;

    
    double T = 0.5; //période de l'onde
    value = 1 + sin(2*M_PI*mytime/T);
    
}