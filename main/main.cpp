#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"
#include "structures.h"

int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    std::vector<std::string> modelNames; // string that contains the name of the models.
    Element mainElement; // The main elements of the mesh.
    Element frontierElement; // The frontier elements of the mesh.

    View mainView; // View of the results.

    gmsh::initialize(argc, argv); // Initialization of gmsh library.
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::open(argv[1]);                          // reads the msh file

    gmsh::model::list(modelNames);

    mainView.name = "MainView";
    mainView.tag = gmsh::view::add(mainView.name);
    mainView.dataType = "NodeData";
    mainView.modelName = modelNames[0];

    std::cout << mainView.modelName << std::endl;

    meshLoader(mainElement, frontierElement); // Initialization of all quantities required.

    solver(mainElement, frontierElement, mainView); // Solving of the PDE with DG-FEM.


    gmsh::finalize(); // Closes gmsh
    return 0;
}
