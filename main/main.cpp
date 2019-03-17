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

    Initialization(mainElements, 2, 3);
    frontierCreation(mainElements, frontierElement, 2, sortedNodes);

    // Loop over all entities.
    

    gmsh::finalize(); // Closes gmsh
    return 0;
}
