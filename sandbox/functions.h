#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Sorting the nodes of the list of the edges to remove duplicate
void mysorting(std::vector<int> & nodes);

// Return the normals at each edges
std::vector<double> my2Dnormal(std::vector<int> nodes);

// Compute the integration of Gauss
void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & basisFunctions,
 const std::vector<double> & determinants, std::vector<double> & massMatrix, const int numElements, const int numGaussPoints);

#endif