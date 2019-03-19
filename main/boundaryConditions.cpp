#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void boundaryConditions(std::vector<double> & coord,double time, double & value){
    
    double L = 0.25; //longueur d'onde
    value=sin(2*M_PI*time); // vague comme quand une goute tombe dans l'eau
}