#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"
#include "structures.h"

// Function that computes the normal to each element of the frontier.

void normals(Element & frontierElement){

    std::size_t i, j;

    // 2D mesh case. The normal is simply the vectorial product of the gradient of the shape
    // function (which is perpendicular to the element edge) and the normal to the plane
    // of the 2D mesh (that is a normal vector along z).
    if(frontierElement.dim == 1)
        for(i = 0; i < frontierElement.elementTag.size(); ++i){
            for(j = 0; j < frontierElement.gaussType; ++j)

            int gradIndex = 3 * frontierElement.neighbours[i].first * \
                           frontierElement.numComponentShapeFunctionsGrad *\
                           frontierElement.gaussType;

            double norm = sqrt(frontierElement.shapeFunctionsGrad[gradIndex + 1] *\
                               frontierElement.shapeFunctionsGrad[gradIndex + 1] +\
                               frontierElement.shapeFunctionsGradParam[gradIndex] *
                               frontierElement.shapeFunctionsGradParam[gradIndex]);

            frontierElement.normals.push_back(frontierElement.shapeFunctionsGradParam[gradIndex + 1]/norm);
            frontierElement.normals.push_back(-frontierElement.shapeFunctionsGradParam[gradIndex]/norm);
            frontierElement.normals.push_back(0.);

        }

    // 3D case mesh. gmsh directly implements the normal to an element face.
    else if(frontierElement.dim == 2)
        for(i = 0; i < frontierElement.neighbours.size(); ++i){

            int gaussIndex = 3 * frontierElement.neighbours[i].first * frontierElement.gaussType;

            std::vector<double> normalTmp;
            std::vector<double> gaussTmp = {frontierElement.gaussPointsParam[gaussIndex],\
                                            frontierElement.gaussPointsParam[gaussIndex + 1]};

            gmsh::model::getNormal(frontierElement.elementTag[i], gaussTmp, normalTmp);

            for(j = 0; j < normalTmp.size(); ++j) frontierElement.normals.push_back(normalTmp[j]);

        }   

}
