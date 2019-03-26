#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void boundaryConditions(std::vector<double> & coord,double mytime, double & value){
    
    //double T = 2; //longueur d'onde
    //value=sin(2*M_PI*mytime/T);
    if(coord[0] == 0)
        value = 0;
    else
        value = 0;
}