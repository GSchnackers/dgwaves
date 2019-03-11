#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

void (std::<double> coord, double & value){
    
    L=0.25; //longueur d'onde
    value=sin(sqrt(coord[1]*coord[1]+coord[2]*coord[2])*2*pi/L)/sqrt(coord[1]*coord[1]+coord[2]*coord[2]); // vague comme quand une goute tombe dans l'eau
}