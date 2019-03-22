#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void initialCondition(std::vector<double> & coord, double & value){
    
    double L=0.25; //longueur d'onde

    value=0;

    //value=sin(sqrt(coord[1]*coord[1]+coord[2]*coord[2])*2*M_PI/L)/sqrt(coord[1]*coord[1]+coord[2]*coord[2]); // vague comme quand une goute tombe dans l'eau

}