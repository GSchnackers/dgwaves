#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

int main(int argc, char **argv)
{

    struct Element element2D; // Structure meant to contain the properties of the elements of dimension 2.

    gmsh::vectorpair entities2D; // Vector pair used to contain the dimension and tags of entities.

    std::size_t i, j; // Index variables

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    gmsh::initialize(argc, argv); // Initialization of gmsh library.
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::open(argv[1]);                            // reads the msh file

    elementInitialization(element2D); // Initialization of the properties of the 2D elements of the same type.

    gmsh::model::getEntities(entities2D, 2); // Get the tag of all geometrical entities of dimension 2.

    std::vector<struct Entity> geometry(entities2D.size()); // Table of structures containing the different entities.

    // Loop over all entities.
    for(i = 0; i < entities2D.size(); ++i){

        geometry[i].entityTag = entities2D[i].second;
        gmsh::model::mesh::getElementsByType(element2D.elementTypes[0], geometry[i].elementTag, \
                                             geometry[i].nodeTags, geometry[i].entityTag);
        gmsh::model::mesh::getElementEdgeNodes(element2D.elementTypes[0], \
                                               geometry[i].elementEdgeNodes, geometry[i].entityTag);
        edges(element2D.numNodes, geometry[i]);

    }

    gmsh::finalize(); // Closes gmsh
    return 0;
}
