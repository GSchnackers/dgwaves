// Implementation of Runge Kutta of order 4

// Input: - vector u at time step t
//        - time step
//        - derivative of u

// Output: -  vector u at time step t+1

#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"

void RungeKutta4(std::vector<int> & nodeTags2D, std::vector<int> & nodeTags2DPlusBC, std::vector<double> & nodeCoord,\
            std::vector<double> & nodeCoordParam, double mytime, double value, std::vector<double> & matrixF,\
            std::vector<double> & uPlusBC, std::vector<double> & u, std::vector<double> & matrixS, double timeStep,\
            std::vector<int> & tagElement1DSorted, std::vector<int> & upwind, std::vector<int> & neighbours1D,\
            std::vector<int> & indicesNei1, std::vector<int> & indicesNei2, std::vector<double> & matrixM_Inverted,\
            std::vector<double> & dudt, std::vector<int> & elementTags2D, int numNodes2D, int NumNodesSide){
    
    std::vector<double> u1(nodeTags2D.size());
    std::vector<double> u2(nodeTags2D.size());
    std::vector<double> u3(nodeTags2D.size());

    std::vector<double> dudt1(u.size());
    std::vector<double> dudt2(u.size());
    std::vector<double> dudt3(u.size());

    for(std::size_t i = 0; i < u.size(); i++){
        u1[i] = u[i];
        u2[i] = u[i];
        u3[i] = u[i];
    }
    slope(nodeTags2D, nodeTags2DPlusBC, nodeCoord, nodeCoordParam, mytime, value, matrixF, uPlusBC, u, matrixS,\
                tagElement1DSorted, upwind, neighbours1D, indicesNei1, indicesNei2, matrixM_Inverted, dudt,\
                elementTags2D, numNodes2D, NumNodesSide); // slope1
    Forward_Euler_method(u1, timeStep/2, dudt); // u_(n+1) = u_n + h/2 * slope1
    slope(nodeTags2D, nodeTags2DPlusBC, nodeCoord, nodeCoordParam, mytime, value, matrixF, uPlusBC, u1, matrixS,\
                tagElement1DSorted, upwind, neighbours1D, indicesNei1, indicesNei2, matrixM_Inverted, dudt1,\
                elementTags2D, numNodes2D, NumNodesSide); // slope2
    Forward_Euler_method(u2, timeStep/2, dudt1); // u_(n+1) = u_n + h/2 * slope2
    slope(nodeTags2D, nodeTags2DPlusBC, nodeCoord, nodeCoordParam, mytime, value, matrixF, uPlusBC, u2, matrixS,\
                tagElement1DSorted, upwind, neighbours1D, indicesNei1, indicesNei2, matrixM_Inverted, dudt2,\
                elementTags2D, numNodes2D, NumNodesSide); // slope3
    Forward_Euler_method(u3, timeStep, dudt2); // u_(n+1) = u_n + h * slope3
    slope(nodeTags2D, nodeTags2DPlusBC, nodeCoord, nodeCoordParam, mytime, value, matrixF, uPlusBC, u3, matrixS,\
                tagElement1DSorted, upwind, neighbours1D, indicesNei1, indicesNei2, matrixM_Inverted, dudt3,\
                elementTags2D, numNodes2D, NumNodesSide); // slope4
    
    for(std::size_t i = 0; i < dudt1.size(); i++)
        dudt[i] += 2*dudt1[i] + 2*dudt2[i] + dudt3[i];
    
    
    Forward_Euler_method(u, timeStep/6, dudt);
}