#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

void meshLoader(Element & mainElements, Element & frontierElement){

    std::size_t i, j;
    std::vector<int> sortedNodes; // Vector of nodes that serves as a basis for the creation of all 1D elements.

    // Initialization of the elements of the mesh.
    gmsh::logger::write("Initializing the main elements of the mesh...");
    Initialization(mainElements, 2, "Gauss4");
    std::cout << "Done." << std::endl;

    // Sorting of the nodes at the frontier of each element.
    gmsh::logger::write("Collection of the information for creating the frontier elements...");
    sortingNeighbouring(mainElements, frontierElement, sortedNodes);
    std::cout << "Done." << std::endl;

    // Creation of the frontier elements on the basis of the vector of sorted nodes.
    gmsh::logger::write("Creation of the frontier elements...");
    frontierCreation(mainElements, frontierElement, 2, sortedNodes);
    std::cout << "Done." << std::endl;

    // Initialization of the element representing the frontiers.
    gmsh::logger::write("Initialization of the frontier elements...");
    Initialization(frontierElement, 1, "Gauss3", true);
    std::cout << "Done." << std::endl;

    // getting the normals of the edges.
    gmsh::logger::write("Computation of the normals at the frontier elements...");
    normals(frontierElement);
    std::cout << "Done." << std::endl;

    // Correspondance computation between the nodes of each frontier element and its index in the general indexations.
    // This function links the nodes of the frontier elements with their indices in the global numerotation.
    gmsh::logger::write("Correspondace computation...");
    correspondance(mainElements, frontierElement);
    std::cout << "Done." << std::endl;

    // Matrix M of the elements loading.
    gmsh::logger::write("Computation of the mass matrix of each element...");
    matrixMaker(mainElements, "M");
    std::cout << "Done." << std::endl;

    // Mass matrix inversion.
    gmsh::logger::write("Computation of the inverse of the mass matrix of each element...");
    for(i = 0; i < mainElements.massMatrix.size(); i += mainElements.numNodes * mainElements.numNodes)
    {
        std::vector<double> tmp(mainElements.numNodes * mainElements.numNodes);
        std::vector<double> inverseTmp(mainElements.numNodes * mainElements.numNodes);

        for(j = 0; j < mainElements.numNodes * mainElements.numNodes; ++j)
            tmp[j] = mainElements.massMatrix[i + j];

        invert(tmp, inverseTmp);

        mainElements.massMatrixInverse.resize(mainElements.massMatrix.size());

        for(j = 0; j < mainElements.numNodes * mainElements.numNodes; ++j)
            mainElements.massMatrixInverse[i + j] = inverseTmp[j];

    }
    std::cout << "Done." << std::endl;

    // Matrix S of the elements loading.
    gmsh::logger::write("Computation of the stiffness matrix of each element...");
    matrixMaker(mainElements, "SX");
    if(mainElements.dim >= 2) matrixMaker(mainElements, "SY");
    else mainElements.stiffnessMatrixY.resize(mainElements.stiffnessMatrixX.size());
    if(mainElements.dim == 3) matrixMaker(mainElements, "SZ");
    else mainElements.stiffnessMatrixZ.resize(mainElements.stiffnessMatrixX.size());
    std::cout << "Done." << std::endl;

    // Setting of the boundary types.
    gmsh::logger::write("Setting the boundary condition type...");
    setBoundaryConditions(frontierElement);
    std::cout << "Done." << std::endl;

}