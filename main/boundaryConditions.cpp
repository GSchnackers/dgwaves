#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void boundaryConditions(std::vector<double> & coord,double time, double & value){
    
    double T = 2; //longueur d'onde
    value=sin(2*M_PI*time/T);
}