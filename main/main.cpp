#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

int main(int argc, char **argv)
{

    std::size_t i, j;

    std::vector<int> sortedNodes; // Vector of nodes that serves as a basis for the creation of all 1D elements.
    std::vector<int> normals; // Vector containing the normal to each frontier element.
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
    Initialization(mainElements, 2, 3);

    // Sorting of the nodes at the frontier of each element.
    sortingNeighbouring(mainElements, frontierElement, sortedNodes);

    // Creation of the frontier elements on the basis of the vector of sorted nodes.
    frontierCreation(mainElements, frontierElement, 2, sortedNodes);

    // Initialization of the element representing the frontiers.
    Initialization(frontierElement, 1, 3);


    gmsh::finalize(); // Closes gmsh
    return 0;
}
