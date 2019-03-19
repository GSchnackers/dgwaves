#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <Eigen/Dense>
#include <cmath>
#include "functions.h"

// Careful, matrix contains only 1 matrix!!!

void invert(std::vector<double> matrix, std::vector<double> & inverse){

    std::size_t i, j;
    int lineSize = std::sqrt(matrix.size());
    Eigen::MatrixXd tmp(lineSize, lineSize);

    for(i = 0; i < matrix.size(); i += lineSize)
        for(j = 0; j < lineSize; ++j)
            tmp(i,j) = matrix[lineSize*i + j];

    tmp = tmp.inverse();

    for(i = 0; i < matrix.size(); i += lineSize)
        for(j = 0; j < lineSize; ++j)
            inverse.push_back(tmp(i,j));

}