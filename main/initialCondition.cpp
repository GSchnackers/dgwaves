#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void initialCondition(std::vector<double> & coord, double & value){
    
    value = 0;

    /*
    double xCentre=0.5;
    double yCentre=0.5;
    double L=0.8; //longueur d'onde
    double R = sqrt((xCentre-coord[0])*(xCentre-coord[0])+(yCentre-coord[1])*(yCentre-coord[1]));
    value = 1 + cos(2*M_PI*R/L)/(1+R); // vague comme quand une goute tombe dans l'eau
    */
}