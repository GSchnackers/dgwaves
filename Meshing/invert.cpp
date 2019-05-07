#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

void invert(std::vector<double> matrix, std::vector<double> & inverse){

    std::size_t i, j; // loop variables.

    // Initialization of the inverse size
    inverse.resize(matrix.size());

    // Computation of the number of rows (= number of columns as the matrix is square).
    int lineSize = std::sqrt(matrix.size());

    // Matrix that is going to be inverted.
    Eigen::MatrixXd tmp(lineSize, lineSize);

    // Transfusion from the vector to the matrx type.
    for(i = 0; i < lineSize; ++i)
        for(j = 0; j < lineSize; ++j)
            tmp(i,j) = matrix[i*lineSize + j];


    // Inversion
    tmp = tmp.inverse();

    // Transfusion from the matrix to the inverted vector.
    for(i = 0; i < lineSize; ++i)
        for(j = 0; j < lineSize; ++j)
            inverse[lineSize * i + j] = tmp.coeff(i,j);
            
}