#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Sorting the nodes of the list of the edges to remove duplicate
void mysorting(std::vector<int> & nodes);

// Return the normals at each edges
std::vector<double> my2Dnormal(std::vector<int> nodes);

#endif