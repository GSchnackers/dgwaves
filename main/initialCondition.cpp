#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void initialCondition(std::vector<double> & coord, double & value){
    
    double L=0.25; //longueur d'onde
    //double R = sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
    value=1.0;
    //value = sin(2*M_PI*R/L)/R; // vague comme quand une goute tombe dans l'eau

}