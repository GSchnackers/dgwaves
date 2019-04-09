#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void initialCondition(std::vector<double> & coord, double & value){
    
    value = 1;

    /*
    double L = 2; //longueur d'onde
    double k = 0.5;
    double l = 0.5;
    double R = sqrt((coord[0]-k)*(coord[0]-k)+(coord[1]-l)*(coord[1]-l));
    value = 1 + cos(R/L)/(0.1+R); // vague comme quand une goute tombe dans l'eau
    */
}