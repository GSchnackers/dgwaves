#include <cstdio>
#include <iostream>
#include <vector>
#include "structures.hpp"

// Function that computes the normal to each element of the frontier. 
// The normal is always outward pointing with respect to the first element of neighbours.

void normals(Element & frontierElement, Element & mainElement){

    std::size_t i, j, k;

    // Temporary vector for the normals.
    std::vector<double> tmpNorm(frontierElement.elementTag.size() * frontierElement.numGp * 3);

    // 2D mesh case. The normal is simply the vectorial product of the gradient of one of the shape
    // function (which is perpendicular to the element edge) and the normal to the plane
    // of the 2D mesh (that is a normal vector along z). We take here the first shape function, as it
    // decays in such a way that the vectorial product between the gradient and the z axis is outward pointing.

    for(i = 0; i < frontierElement.elementTag.size(); ++i) // Run through the elements
        for(j = 0; j < frontierElement.numGp; ++j) // run through the gauss points of a given element.
        {
            
            int jacobianIndex = i *  frontierElement.numGp * 9 + j * 9;

            int frontierIndex = i * frontierElement.numGp * 3 + j * 3;

            int mainJacobIndex = frontierElement.neighbours[i].first * mainElement.numGp * 9; 

            double compoX = 0, compoY = 0, compoZ = 0;

            double sign = 1;

            double norm;

            if(frontierElement.dim == 1)
            {
                compoX = frontierElement.jacobiansInverse[jacobianIndex + 1];
                compoY = frontierElement.jacobiansInverse[jacobianIndex + 4];

                sign = mainElement.jacobiansInverse[mainJacobIndex + 8];
            }

            else if(frontierElement.dim == 2)
            {
                compoX = frontierElement.jacobiansInverse[jacobianIndex + 2];
                compoY = frontierElement.jacobiansInverse[jacobianIndex + 5];
                compoZ = frontierElement.jacobiansInverse[jacobianIndex + 8];
            }

            norm = sqrt(compoX * compoX + compoY * compoY + compoZ * compoZ);
            
            tmpNorm[frontierIndex] = sign * compoX/norm;
            tmpNorm[frontierIndex + 1] = sign * compoY/norm;
            tmpNorm[frontierIndex + 2] = sign * compoZ/norm;

        }


    frontierElement.normals = tmpNorm; 

}
