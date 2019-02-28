// my2Dnormal is a function that computes the normal to each edge in 2D.
// Input: nodes, the vector containing all nodes (no doublons)
// Output: a vector containing the components of the normal to the edge in the
// format (n1x, n1y, n2x, n2y, n3x, n3y, ...).
// This method is inneficient but works as a first approximation for a preliminary solver.

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include "functions.h"

std::vector<double> my2Dnormal(std::vector<int> nodes)
{

    // Declaration of the vector contaning the normals to the edges.
    std::vector<double> normal2D(nodes.size());

    for(std::size_t i = 0; i < nodes.size(); i += 2)
    {
        std::vector<double> nodeCoord1, nodeCoord2, nodeCoordParam1, nodeCoordParam2;

        gmsh::model::mesh::getNode(nodes[i], nodeCoord1, nodeCoordParam1);
        gmsh::model::mesh::getNode(nodes[i + 1], nodeCoord2, nodeCoordParam2);

        // Computation of the normal. n = (-y , x)/(x^2+y^2)^(1/2)
        
        normal2D[i] = nodeCoord1[1] - nodeCoord2[1]; // -y
        normal2D[i + 1] = nodeCoord2[0] - nodeCoord1[0]; // x

        double norm = sqrt(normal2D[i] * normal2D[i] + normal2D[i + 1] * normal2D[i + 1]); // x^2 + y^2

        // Final norm.
        normal2D[i] /= norm;
        normal2D[i + 1] /= norm;

    }

    return normal2D;

}

/* getNode

    Get the coordinates and the parametric coordinates (if any) of the node with tag tag. This is a sometimes useful but inefficient way of accessing nodes, as it relies on a cache stored in the model. For large meshes all the nodes in the model should be numbered in a continuous sequence of tags from 1 to N to maintain reasonnable performance (in this case the internal cache is based on a vector; otherwise it uses a map).

    Input:

        nodeTag 
    Output:

        coord, parametricCoord 
    Return:

        - 
 */