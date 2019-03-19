#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Forward Euler method
void Forward_Euler_method(std::vector<double> & u, const double timestep, const std::vector<double> & dudt);

// Compute the integration of Gauss
void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & functions,
 const std::vector<double> & determinants, std::vector<double> & matrix,
  const int numElements, const int numGaussPoints, const int numNodes);

// Computes the gradient of a lagrangian shape function at one point.
void gradient(const double f,const  std::vector<double> pointsPositions,\
                             const std::vector<double> nodesPositions, std::vector<double> & grad);

// Allow to compute the neighbours of each edges.
void neighbours(const std::vector<int> nodeTags, const int nodeNumber,const std::vector<int> elementTags,\
                const std::vector<int> nodes, std::vector<int> & neighbourhood);

// Return the normals at each edges
void normal(const std::vector<int> nodes, std::vector<double> & normal2D);

// Sorting the nodes of the list of the edges to remove duplicate
void sorting(std::vector<int> & nodes);

// Vector of flux
void vectorF(std::vector<double> & vectorF, const double a_x, const double a_y,\
                    const std::vector<double> u, const std::vector<double> sfg);

// Initial condition
void initialCondition(std::vector<double> & coord, double & value);

// Boundary conditions
void boundaryConditions(std::vector<double> & coord,double time, double & value);

// Inversion of a square dense matrix.
void invert(std::vector<double> matrix, std::vector<double> & inverse);

#endif