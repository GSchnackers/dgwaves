#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "meshing.hpp"
#include "structures.hpp"

void meshLoader(Element & mainElements, Element & frontierElement, std::string & gaussType, \
                PhysicalGroups & physicalGroups, int mainDim){

    std::size_t i, j;
    std::vector<int> sortedNodes; // Vector of nodes that serves as a basis for the creation of all 1D elements.

    int frontierDim = mainDim - 1;

    // Initialization of the elements of the mesh.
    gmsh::logger::write("Initializing the main elements of the mesh...");
    Initialization(mainElements, mainDim, gaussType);
    gmsh::logger::write("Done.");

    // Sorting of the nodes at the frontier of each element.
    gmsh::logger::write("Collection of the information for creating the frontier elements...");
    sortingNeighbouring(mainElements, frontierElement, sortedNodes);
    gmsh::logger::write("Done.");

    // Creation of the frontier elements on the basis of the vector of sorted nodes.
    gmsh::logger::write("Creation of the frontier elements...");
    frontierCreation(mainElements, frontierElement, mainDim, sortedNodes);
    gmsh::logger::write("Done.");

    // Initialization of the element representing the frontiers.
    gmsh::logger::write("Initialization of the frontier elements...");
    Initialization(frontierElement, frontierDim, gaussType, true);
    gmsh::logger::write("Done.");

    // getting the normals of the edges.
    gmsh::logger::write("Computation of the normals at the frontier elements...");
    normals(frontierElement, mainElements);
    gmsh::logger::write("Done.");

    /*
    for(i = 0; i < frontierElement.normals.size(); ++i)
        std::cout << frontierElement.normals[i] << " " << mainElements.elementTag[frontierElement.neighbours[i/(3 * frontierElement.numGp)].first] << std::endl;
    */

    // Correspondance computation between the nodes of each frontier element and its index in the general indexations.
    // This function links the nodes of the frontier elements with their indices in the global numerotation.
    gmsh::logger::write("Correspondace computation...");
    correspondance(mainElements, frontierElement);
    gmsh::logger::write("Done.");

    // Matrix M of the elements loading.
    gmsh::logger::write("Computation of the mass matrix of each element...");
    matrixMaker(mainElements, "M");
    gmsh::logger::write("Done.");

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
    gmsh::logger::write("Done.");

    
    // Matrix S of the elements loading.
    gmsh::logger::write("Computation of the stiffness matrix of each element...");
    matrixMaker(mainElements, "SX");
    if(mainElements.dim >= 2) matrixMaker(mainElements, "SY");
    else mainElements.stiffnessMatrixY.resize(mainElements.stiffnessMatrixX.size(), 0);
    if(mainElements.dim == 3) matrixMaker(mainElements, "SZ");
    else mainElements.stiffnessMatrixZ.resize(mainElements.stiffnessMatrixX.size(), 0);
    gmsh::logger::write("Done.");

    // Loading of all physical Groups.
    gmsh::logger::write("Computation of the physical groups...");

    gmsh::model::getPhysicalGroups(physicalGroups.dimTags);
    physicalGroups.entityTags.resize(physicalGroups.dimTags.size());
    physicalGroups.name.resize(physicalGroups.dimTags.size());
    physicalGroups.elemType.resize(physicalGroups.dimTags.size());
    for(i = 0; i < physicalGroups.dimTags.size(); ++i)
    {
        std::vector<std::vector<int>> bin1, bin2;
        gmsh::model::getPhysicalName(physicalGroups.dimTags[i].first, \
                                     physicalGroups.dimTags[i].second,\
                                     physicalGroups.name[i]);

        gmsh::model::getEntitiesForPhysicalGroup(physicalGroups.dimTags[i].first, \
                                                 physicalGroups.dimTags[i].second, \
                                                 physicalGroups.entityTags[i]);

        physicalGroups.elemType[i].resize(physicalGroups.entityTags[i].size());

        for(j = 0; j < physicalGroups.entityTags[i].size(); ++j)
            gmsh::model::mesh::getElements(physicalGroups.elemType[i][j], bin1 , bin2, \
                                           physicalGroups.dimTags[i].first, physicalGroups.entityTags[i][j]);
    }
    
    gmsh::logger::write("Done.");

}