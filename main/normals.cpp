#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"
#include "structures.h"

// Function that computes the normal to each element of the frontier. 
// The normal is always outward pointing with respect to the first element of neighbours.

void normals(Element & frontierElement){

    std::size_t i, j, k;

    // Temporary vector for the normals.
    std::vector<double> tmpNorm(frontierElement.elementTag.size() * frontierElement.numGp * 3);

    // 2D mesh case. The normal is simply the vectorial product of the gradient of one of the shape
    // function (which is perpendicular to the element edge) and the normal to the plane
    // of the 2D mesh (that is a normal vector along z). We take here the first shape function, as it
    // decays in such a way that the vectorial product between the gradient and the z axis is outward pointing.

    if(frontierElement.dim == 1)
        for(i = 0; i < frontierElement.elementTag.size(); ++i) // Run through the elements
            for(j = 0; j < frontierElement.numGp; ++j) // run through the gauss points of a given element.
            {

                int jacobianIndex = i *  frontierElement.numGp * 9 + j * 9;

                int frontierIndex = i * frontierElement.numGp * 3 + j * 3;

                double norm = sqrt(frontierElement.jacobiantInverseTranspose[jacobianIndex + 3] *\
                                frontierElement.jacobiantInverseTranspose[jacobianIndex + 3] +\
                                frontierElement.jacobiantInverseTranspose[jacobianIndex] *
                                frontierElement.jacobiantInverseTranspose[jacobianIndex]); // Normalization.

                int jacobIndex = frontierElement.neighbours[frontierIndex].first;

                tmpNorm[frontierIndex] = frontierElement.jacobiantInverseTranspose[jacobianIndex + 3]/norm;
                tmpNorm[frontierIndex + 1] = -frontierElement.jacobiantInverseTranspose[jacobianIndex ]/norm;
                tmpNorm[frontierIndex + 2] = 0.;

            }

    // 3D case mesh. gmsh directly implements the normal to an element face.
    else if(frontierElement.dim == 2)
        for(i = 0; i < frontierElement.neighbours.size(); ++i){

            int gaussIndex = 3 * frontierElement.neighbours[i].first * frontierElement.numGp;

            std::vector<double> normalTmp;
            std::vector<double> gaussTmp = {frontierElement.gaussPointsParam[gaussIndex],\
                                            frontierElement.gaussPointsParam[gaussIndex + 1]};

            gmsh::model::getNormal(frontierElement.elementTag[i], gaussTmp, normalTmp);

            for(j = 0; j < normalTmp.size(); ++j) tmpNorm[3*i + j] = normalTmp[j];

        } 


    frontierElement.normals = tmpNorm; 

}
