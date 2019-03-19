#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <Eigen/Dense>
#include <cmath>
#include "functions.h"

// Carefus, matrix contains only 1 matrix!!!

void invert(std::vector<double> matrix, std::vector<double> & invert){

    std::size_t i, j;
    Eigen::MatrixXd tmp(int(std::sqrt(matrix.size())),int(std::sqrt(matrix.size())));
    int lineSize = std::sqrt(matrix.size());

    for(i = 0; i < matrix.size(); ++i){

        //tmp();

    }

}