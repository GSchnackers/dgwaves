#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

int main(int argc, char **argv)
{

    std::size_t i, j;

    std::vector<int> sortedNodes; // Vector of nodes that serves as a basis for the creation of all 1D elements.
    struct Element mainElements; // The main elements of the mesh.
    struct Element frontierElement; // The frontier elements of the mesh.

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    gmsh::initialize(argc, argv); // Initialization of gmsh library.
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::open(argv[1]);                            // reads the msh file

    // Initialization of the elements of the mesh.
    gmsh::logger::write("Initializing the main elements of the mesh...");
    Initialization(mainElements, 2, "Gauss3");
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

    for(i = 0; i< frontierElement.normals.size(); ++i)
        std::cout<< frontierElement.normals[i] << std::endl;
    

    gmsh::finalize(); // Closes gmsh
    return 0;
}
