#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Sorting the nodes of the list of the edges to remove duplicate
void mysorting(std::vector<int> & nodes);

// Return the normals at each edges
void normal(const std::vector<int> nodes, std::vector<double> & normal2D);

// Compute the integration of Gauss
void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & functions,
 const std::vector<double> & determinants, std::vector<double> & matrix, const int numElements, const int numGaussPoints);

// Computes the gradient of a lagrangian shape function at one point.

void gradient(const double f,const  std::vector<double> pointsPositions,\
                             const std::vector<double> nodesPositions, std::vector<double> & grad);

#endif