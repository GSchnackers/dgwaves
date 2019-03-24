#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

// Matrix maker builds any stiffness or mass matrix that exists provided the element and the matrix type.
// The matrix builded are assumed to be symmetric.

void matrixMaker(Element & element, std::string matrixType)
{

    std::size_t i, j, k, l; // Indiex variable.

    // Matrix containing all values for all elements and all gauss points.
    std::vector<double> tmp1 = element.shapeFunctionsParam, tmp2;
    std::vector<double> matrixTmp(element.numNodes * element.numNodes * element.elementTag.size(), 0);

    int compo = -1; // Keeps in memory the component of the gradient to use in the building of x. It is -1 for M, 0, 1 and 2 for SX SY and SZ respectively.

    if(!matrixType.compare("M"))
        tmp2 = tmp1;

    else if(!matrixType.compare("SX") || !matrixType.compare("SY") || !matrixType.compare("SZ"))
    {

        if(matrixType.find("X") != std::string::npos) compo = 0;
        else if(matrixType.find("Y") != std::string::npos) compo = 1;
        else compo = 2;

        tmp2.resize(tmp1.size());

        for(i = compo; i < element.shapeFunctionsGradParam.size(); i += 3) 
            tmp2[i/3] = element.shapeFunctionsGradParam[i];
            
    }

    else
    {
        gmsh::logger::write("There matrix type is not recognized.", "error");
        exit(-1);
    }

    for(i = 0; i < element.elementTag.size(); ++i)
        for(j = 0; j < element.numNodes; ++j)
            for(k = j; k < element.numNodes; ++k)
            {
                int indexMatrix = i * element.numNodes * element.numNodes + j * element.numNodes + k;
                int indexMatrixTranspose = i * element.numNodes * element.numNodes + k * element.numNodes + j;

                for(l = 0; l < element.numGp; ++l)
                {
                    int index1 = j * element.numGp + l;
                    int index2 = k * element.numGp + l;
                    int indexJacob = i * element.numGp + l;
                    int indexGPoint1 = j * element.numGp + 3;
                    int indexGPoint2 = k * element.numGp + 3;

                    matrixTmp[indexMatrix] += tmp1[index1] * tmp2[index2] * \
                                              element.jacobiansDet[indexJacob] * \
                                              element.gaussPointsParam[indexGPoint1] * \
                                              element.gaussPointsParam[indexGPoint2];
                    
                }

                matrixTmp[indexMatrixTranspose] = matrixTmp[indexMatrix];

            }

    switch (compo)
    {
        case -1:
            element.massMatrix = matrixTmp;
            break;

        case 0:
            element.stiffnessMatrixX = matrixTmp;
            break;
        
        case 1:
            element.stiffnessMatrixY = matrixTmp;
            break;

        case 2:
            element.stiffnessMatrixZ = matrixTmp;
            break;
    
    }
    
}