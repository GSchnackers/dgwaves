// Implementation of the matrix F

/* Input: 
- coef a_x and a_y
- values of u at the nodes
- values of the basis function at gauss points of the edge (bfg = basis functions at gauss point) (size = N_basis_function * N_gausspoint)
- edgeNodes
*/   

// Output: -  matrix F (size N_elements * N_basis_function) attention, il faut l'initialiser à 0 après chaque pas de temps

#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"
#include <gmsh.h>

void matrixF(std::vector<double> & matrixF, std::vector<double> coefTransport){




}