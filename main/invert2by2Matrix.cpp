#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

/*
INPUT :
- the matrix to be inverted : Matrix = [comp_{11}, comp_{12},...,comp_{22}]

OUTPUT :
- the inverse of the matrix : MatrixInverted = [comp_{11}, comp_{12},...,comp_{22}]
*/



invert2by2Matrix(std::vector<double> Matrix, std::vector<double> & MatrixInverted){

double det = Matrix[0]*Matrix[3]-Matrix[1]*Matrix[2];

MatrixInverted[0] = Matrix[3]/det;
MatrixInverted[1] = - Matrix[1]/det;
MatrixInverted[2] = - Matrix[2]/det;
MatrixInverted[3] = Matrix[0]/det;


return 
}