#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

int main(int argc, char **argv)
{

    std::size_t i, j;
    std::vector<struct Entity> geometry; // Structure containing all info about the model.

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    gmsh::initialize(argc, argv); // Initialization of gmsh library.
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::open(argv[1]);                            // reads the msh file
    Initialization(geometry);

    for(i = 0; i < geometry.size(); ++i){

        std::cout << geometry[i].entityTag << std::endl;
        std::cout << std::endl;
        for(j = 0; j < geometry[i].nodeTags.size(); ++j) std::cout << geometry[i].nodeTags[j] << std::endl;
        std::cout << std::endl;
        for(j = 0; j < geometry[i].elementTag2D.size(); ++j) std::cout << geometry[i].elementTag2D[j] << std::endl;
        std::cout << std::endl;
        for(j = 0; j < geometry[i].elementEdgeNode2D.size(); ++j) std::cout << geometry[i].elementEdgeNode2D[j] << std::endl;
        std::cout << std::endl;

    }

    // Loop over all entities.
    

    gmsh::finalize(); // Closes gmsh
    return 0;
}
