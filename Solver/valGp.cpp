/* This functions computes the approximated function (Galerkin) at the Gauss points of the elements
   described by "element" in parametric coordinates. It stores the results in the "result" vector in the form
   [e1G1, e1G2, ... , e2G1, e2G2, ...]. */

#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "structures.hpp"

double valGpBound(int i, int gpIndex, double t, const Quantity & u, const Element & frontierElement, \
                  const Properties & matProp){

    int paramIdx = 9 * i;

    if(frontierElement.neighbours[i].second == OPENING)
        return u.gp[gpIndex].first;

    else if(frontierElement.neighbours[i].second == SINUS_E)
    {

        if(gpIndex % 6 == 0) 
            return frontierElement.bcParam[paramIdx] * sin(frontierElement.bcParam[paramIdx + 1] \
                    * M_PI * t + frontierElement.bcParam[paramIdx + 2]);
                    
        else if(gpIndex % 6 == 1)                           
            return frontierElement.bcParam[paramIdx + 3] * sin(frontierElement.bcParam[paramIdx + 4] \
                    * M_PI * t + frontierElement.bcParam[paramIdx + 5]);

        else if(gpIndex % 6 == 2)                                
            return frontierElement.bcParam[paramIdx + 6] * sin(frontierElement.bcParam[paramIdx + 7] \
                    * M_PI * t + frontierElement.bcParam[paramIdx + 8]); 
                                            
        else return u.gp[gpIndex].first;

    }

    else if(frontierElement.neighbours[i].second == SINUS_H)
    {
        
        if(gpIndex % 6 == 3) 
            return frontierElement.bcParam[paramIdx] * sin(frontierElement.bcParam[paramIdx + 1] \
                    * M_PI * t + frontierElement.bcParam[paramIdx + 2]);
                    
        else if(gpIndex % 6 == 4)                           
            return frontierElement.bcParam[paramIdx + 3] * sin(frontierElement.bcParam[paramIdx + 4] \
                    * M_PI * t + frontierElement.bcParam[paramIdx + 5]);

        else if(gpIndex % 6 == 5)                                
            return frontierElement.bcParam[paramIdx + 6] * sin(frontierElement.bcParam[paramIdx + 7] \
                    * M_PI * t + frontierElement.bcParam[paramIdx + 8]); 
                                            
        else return u.gp[gpIndex].first;

    }

    else if(frontierElement.neighbours[i].second == PERFECTCOND)
    {
        if((gpIndex % 6) < 3) return 0;
        else return u.gp[gpIndex].first;
    }

    else if(frontierElement.neighbours[i].second == TE2D)
    {
        if(gpIndex % 6 == 2)
            return sin(frontierElement.bcParam[paramIdx] * M_PI * (frontierElement.gaussPoints[gpIndex/6 * 3 + 1] - 0.5)) \
                   * cos(frontierElement.bcParam[paramIdx + 1] * M_PI * t);
                   
        else 
            return u.gp[gpIndex].first;
    }

    else if(frontierElement.neighbours[i].second == TE3D)
    {
        if(gpIndex % 6 == 1)
            return cos(frontierElement.bcParam[paramIdx] * M_PI * (frontierElement.gaussPoints[gpIndex/6 * 3 + 1] - 0.5)) \
                   * sin(frontierElement.bcParam[paramIdx + 1] * M_PI * (frontierElement.gaussPoints[gpIndex/6 * 3 + 2] - 0.125)/0.25) \
                   * cos(frontierElement.bcParam[paramIdx + 2] * M_PI * t);
        
        else if(gpIndex % 6 == 2)
            return sin(frontierElement.bcParam[paramIdx] * M_PI * (frontierElement.gaussPoints[gpIndex/6 * 3 + 1] - 0.5)) \
                   * cos(frontierElement.bcParam[paramIdx + 1] * M_PI * (frontierElement.gaussPoints[gpIndex/6 * 3 + 2] - 0.125)/0.25) \
                   * cos(frontierElement.bcParam[paramIdx + 2] * M_PI * t);
                   
        else 
            return u.gp[gpIndex].first;
    }

    else if(frontierElement.neighbours[i].second == SINE)
    {
        return frontierElement.bcParam[paramIdx] * \
               sin(frontierElement.bcParam[paramIdx + 1] * M_PI * t + frontierElement.bcParam[paramIdx + 2]);
    }

    else
    {
        gmsh::logger::write("Error", "error");
        exit(-1);
    }

}

void valGp(Quantity & u, const Element & mainElement, const Element & frontierElement, int numU, \
           const Properties & matProp, double t){

    std::size_t i, j, k, l;

    std::fill(u.gp.begin(), u.gp.end(), std::make_pair(0,0));

    for(i = 0; i < frontierElement.elementTag.size(); ++i) // Loop over the elements
        for(j = 0; j < frontierElement.numGp; ++j) // Loop over the Gauss Points
            for(k = 0; k < frontierElement.numNodes; ++k) // loop over the nodes
                for(l = 0; l < numU; ++l)
                {
                    int gpIndex    = i * frontierElement.numGp * numU + j * numU + l;
                    int shapeIndex = j * frontierElement.numNodes + k;
    
                    int frontNodeIndex = i * frontierElement.numNodes + k;

                    int mainNodeIndex1 = frontierElement.neighbours[i].first * mainElement.numNodes * numU + \
                                         frontierElement.nodeCorrespondance[frontNodeIndex].first * numU + l;

                    int mainNodeIndex2 = frontierElement.neighbours[i].second * mainElement.numNodes * numU + \
                                         frontierElement.nodeCorrespondance[frontNodeIndex].second * numU + l;
                    

                    u.gp[gpIndex].first += u.node[mainNodeIndex1] * \
                                           frontierElement.shapeFunctionsParam[shapeIndex];
                    

                    if(frontierElement.neighbours[i].second >= 0)
                        u.gp[gpIndex].second += u.node[mainNodeIndex2] * \
                                                frontierElement.shapeFunctionsParam[shapeIndex];

                    else if(k == frontierElement.numNodes - 1)
                        u.gp[gpIndex].second = valGpBound(i, gpIndex, t, u, frontierElement, matProp);
                    
                    
                }
            

}