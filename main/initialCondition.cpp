#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void initialCondition(std::vector<double> & coord, double & value){
    
    value = 1;

    /*
    double L=0.8; //longueur d'onde
    double R = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
    value = 1 + cos(2*M_PI*R/L)/(1+R); // vague comme quand une goute tombe dans l'eau
    */
}