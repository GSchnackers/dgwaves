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
        return -1;
    }

    std::vector<std::string> modelNames; // string that contains the name of the models.
    
    Element mainElement; // The main elements of the mesh.
    Element frontierElement; // The frontier elements of the mesh.

    PhysicalGroups physicalGroups;

    View mainView; // View of the results.

    Simulation simulation; // parameters of the simulation.

    int meshDim; // dimension of the mesh.

    gmsh::initialize(argc, argv); // Initialization of gmsh library.
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::option::setNumber("Mesh.SaveAll", 1);
    gmsh::open(argv[1]); // reads the msh file

    gmsh::logger::write("Simulation parameter loading...");
    readParam(argv[2], simulation);
    gmsh::logger::write("Done.");

    // Gets the dimension of the geometric model, which is the dimension of the mesh.

    meshLoader(mainElement, frontierElement, simulation.gaussType, physicalGroups, gmsh::model::getDimension()); // Initialization of all quantities required.

    gmsh::model::list(modelNames);

    mainView.name = "MainView";
    mainView.tag = gmsh::view::add(mainView.name);
    mainView.dataType = "ElementNodeData";
    
    mainView.modelName = modelNames[0];
    mainView.data.resize(mainElement.elementTag.size(), std::vector<double>(6 * mainElement.numNodes));

    solver(mainElement, frontierElement, physicalGroups, mainView, simulation); // Solving of the PDE with DG-FEM.

    gmsh::finalize(); // Closes gmsh
    return 0;
}
