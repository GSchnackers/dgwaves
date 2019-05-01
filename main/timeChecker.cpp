#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void timeChecker(const Element & mainElement, const Element & frontierElement,\
                 const Quantity & flux, const Quantity & u, const std::vector<double> & SFProd, \
                 const std::vector<double> & fluxVector, const double t){

    std::size_t i, j, k;

    // Boundary Conditions verification.
    std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << " %%%%%%%%%%%%%%%%%%%%% TIME STEP t = " << t << " %%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    std::cout << "BC's verifier at t = " << t << std::endl;
    for(i = 0; i < u.bound.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(6 * mainElement.numNodes)] << " Node: " << mainElement.nodeTags[i/6] << " Value: " << u.bound[i] << " " << std::endl;
    std::cout << std::endl;

    std::cout << "Gauss points verifier at t = " << t << std::endl;
        
    // Value of u at the gauss Point verifier.
    for(i = 0; i < u.gp.size(); ++i)
    {
        std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 6)].first] << " Gauss Point: " << i % (frontierElement.numGp * 6) << " Value: " << u.gp[i].first;
        if(frontierElement.neighbours[i/(frontierElement.numGp * 6)].second >= 0)
            std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 6)].second] << " Gauss Point: " << i % (frontierElement.numGp * 6) << " Value: " << u.gp[i].second << " ";
        else
            std::cout << " / Element: NONE Gauss Point: " << i % (frontierElement.numGp*6) << " Value: " << u.gp[i].second << " ";    
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Value of the nodal flux verifier.
    std::cout << "Physical nodal flux verifier at t = " << t << std::endl;
    for(i = 0; i < flux.node.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(18 * mainElement.numNodes)] << " Node: " << mainElement.nodeTags[i/18] << " Value: " << flux.node[i] << std::endl;
    std::cout << std::endl;

    std::cout << "Physical flux at the Gauss points verifier at t = " << t << std::endl;
    
    // Value of the flux at the gauss point verifyer.
    for(i = 0; i < flux.gp.size(); ++i)
    {
        std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 18)].first] << " Gauss Point: " << (i / 18) % frontierElement.numGp << " Value: " << flux.gp[i].first;
        if(frontierElement.neighbours[i/(18*frontierElement.numGp)].second >= 0)
            std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 18)].second] << " Gauss Point: " << (i / 18) % frontierElement.numGp << " Value: " << flux.gp[i].second << " ";
        else
            std::cout << " / Element: NONE Gauss Point: " << (i / 18) % frontierElement.numGp << " Value: " << flux.gp[i].second << " ";    
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Value of the physical flux verifier.
    /* std::cout << "Physical flux direction at the Gauss points at t = " << t << std::endl;
    for(i = 0; i < flux.direction.size(); ++i)
    {
        std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 18)].first];
        if(frontierElement.neighbours[i/(18*frontierElement.numGp)].second >= 0)
            std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 19)].second] << " Gauss Point: " << (i / 18) % frontierElement.numGp << " Value: " << flux.direction[i] << " ";
        else
            std::cout << " / Element: NONE Gauss Point: " << (i / 18) % frontierElement.numGp << " Value: " << flux.direction[i] << " ";    
        std::cout << std::endl;
    }
    std::cout << std::endl; */ 

    // Value of the numerical flux verifier.
    std::cout << "Numerical flux at the Gauss points verifier at t = " << t << std::endl;
    for(i = 0; i < flux.num.size(); ++i)
    {
        std::cout << "Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 18)].first];
        if(frontierElement.neighbours[i/(18*frontierElement.numGp)].second >= 0)
            std::cout << " / Element: " << mainElement.elementTag[frontierElement.neighbours[i/(frontierElement.numGp * 18)].second] << " Gauss Point: " << (i / 18) % frontierElement.numGp << " Value: " << flux.num[i] << " ";
        else
            std::cout << " / Element: NONE Gauss Point: " << (i / 18) % frontierElement.numGp << " Value: " << flux.num[i] << " ";    
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Stiffness flux product verifier.
    std::cout << "Stiffness/Flux product verifier at t = " << t << std::endl;
    for(i = 0; i < SFProd.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(mainElement.numNodes * 6)] << " Node: " << mainElement.nodeTags[i/6] << " Value: " << SFProd[i] << std::endl;
    std::cout << std::endl;

    // Flux Integration verifier.
    std::cout << "Flux integration verifier at t = " << t << std::endl;
    for(i = 0; i < SFProd.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(mainElement.numNodes * 6)] << " Node: " << mainElement.nodeTags[i/6] << " Value: " << fluxVector[i] << std::endl;
    std::cout << std::endl;

    // Nodal values verifier.
    std::cout << "Nodal value verifier at t = " << t << std::endl;
    for(i = 0; i < u.node.size(); ++i)
        std::cout << "Element: " << mainElement.elementTag[i/(mainElement.numNodes * 6)] << " Node: " << mainElement.nodeTags[i/6] << " Value: " << u.node[i] << std::endl;
    std::cout << std::endl;

}