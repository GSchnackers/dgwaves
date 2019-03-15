#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Initialization of the properties of the element of a certain dim and a certain type.
void Initialization(std::vector<struct Entity> & geometry);

// Sorting the nodes of the list of the edges to remove duplicate and giving the neighbours of each edge.
void edges(const int numNodes, struct Entity & entity);

/*// Return the normals at each edges
void normal(const std::vector<int> nodes, std::vector<double> & normal2D,\
            std::vector<double> & nodeCoords, std::vector<double> & nodeCoordParam);

// Compute the integration of Gauss
void gaussIntegration(const std::vector<double> & integrationPoints, const std::vector<double> & functions,
 const std::vector<double> & determinants, std::vector<double> & matrix,
  const int numElements, const int numGaussPoints, const int numNodes);

// Computes the gradient of a lagrangian shape function at one point.
void gradient(const double f,const  std::vector<double> pointsPositions,\
                             const std::vector<double> nodesPositions, std::vector<double> & grad);

// Allow to compute the neighbours of each edges.

void neighbours(const std::vector<int> nodeTags, const int nodeNumber,\
               const std::vector<int> elementTags, std::vector<int> & nodes); */

#endif