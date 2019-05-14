#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "structures.hpp"

// Matrix maker builds any stiffness or mass matrix that exists provided the element and the matrix type.
// The matrix builded are assumed to be symmetric.

void matrixMaker(Element & element, std::string matrixType)
{

    std::size_t i, j, k, l; // Index variable.

    element.massMatrix.resize(element.elementTag.size() * element.numNodes * element.numNodes);

    if(!matrixType.compare("M"))
    {
        // Computes the resulting matrix.
        for(i = 0; i < element.elementTag.size(); ++i)
            for(j = 0; j < element.numNodes; ++j)
                for(k = 0; k < element.numNodes; ++k)
                {
                    int indexMatrix = i * element.numNodes * element.numNodes + j * element.numNodes + k;

                    for(l = 0; l < element.numGp; ++l)
                        element.massMatrix[indexMatrix] += \
                                            element.shapeFunctionsParam[l * element.numNodes + j] * \
                                            element.shapeFunctionsParam[l * element.numNodes + k] * \
                                            element.jacobiansDet[i * element.numGp + l] * \
                                            element.gaussPointsParam[4 * l + 3];
                    
                }
    } 

    else if(!matrixType.compare("SX") || !matrixType.compare("SY") || !matrixType.compare("SZ"))
    {
        int compo; // Keeps in memory the component of the gradient to use in the building of x. It is 0, 1 and 2 for SX SY and SZ respectively.

        std::vector<double> stiffTmp(element.elementTag.size() * element.numNodes * element.numNodes, 0);
        if(matrixType.find("X") != std::string::npos) compo = 0;
        else if(matrixType.find("Y") != std::string::npos) compo = 1;
        else compo = 2;

        std::vector<double> realGrad(element.elementTag.size() * element.numGp * element.numNodes, 0);

        // Computation of the real gradient.
        for(i = 0; i < element.elementTag.size(); ++i) // loop over the element tags.
            for(j = 0; j < element.numGp; ++j) // loop over the gauss points.
                for(k = 0; k < element.numNodes; ++k) // loop over the nodes.
                {
                    int realIndex = i * element.numGp * element.numNodes + j * element.numNodes + k;

                    for(l = 0; l < 3; ++l) // loop over the lines of the jacobian.
                    {   
                        int jacobIndex = i * element.numGp * 9 + j * 9 + compo * 3 + l;
                        int paramIndex = j * element.numNodes * 3 + k * 3 + l;

                        // Multiplies the corresponding row of the jacobian with the isoparam gradient to obtain the real gradient.
                        realGrad[realIndex] += element.jacobiansInverse[jacobIndex] * \
                                               element.shapeFunctionsGradParam[paramIndex];
                    }
                }

        // Computation of the Stiffness matrix.
        for(i = 0; i < element.elementTag.size(); ++i)
            for(j = 0; j < element.numNodes; ++j)
                for(k = 0; k < element.numNodes; ++k)
                {
                    int stiffIndex =  i * element.numNodes * element.numNodes + j * element.numNodes + k; 

                    for(l = 0; l < element.numGp; ++l)
                    {
                        int realIndex = i * element.numGp * element.numNodes + l * element.numNodes + j;
                        int shapeIndex = l * element.numNodes + k;
                        int detIndex = i * element.numGp + l;

                        stiffTmp[stiffIndex] += realGrad[realIndex] * \
                                                element.shapeFunctionsParam[shapeIndex] * \
                                                element.gaussPointsParam[4 * l + 3] *
                                                element.jacobiansDet[detIndex];
                        
                    }
                }

        switch(compo)
        {
            case 0:
                element.stiffnessMatrixX = stiffTmp;
            break;

            case 1:
                element.stiffnessMatrixY = stiffTmp;
            break;

            case 2:
                element.stiffnessMatrixZ = stiffTmp;
            break;

        }
            
    }

    else
    {
        gmsh::logger::write("The matrix type is not recognized.", "error");
        exit(-1);
    }
    
}