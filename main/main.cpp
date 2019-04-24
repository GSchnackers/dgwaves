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

    View mainView; // View of the results.

    double simTime, timeInc, alpha; // Duration of the simulation and time incrementation.
    int registration, meshDim, solvType, debug; // Type of the solver (0 for euler, 1 for runge-kutta) and registration appear 1/registration steps.
    std::string gaussType; // Guets the integration time.

    gmsh::initialize(argc, argv); // Initialization of gmsh library.
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::option::setNumber("Mesh.SaveAll", 1);
    gmsh::open(argv[1]); // reads the msh file

    gmsh::logger::write("Simulation parameter loading...");
    readParam(argv[2], simTime, timeInc, registration, solvType, gaussType, meshDim, debug, alpha);
    gmsh::logger::write("Done.");

    std::cout << simTime << " " << timeInc << " " << registration << " " << solvType << " " << gaussType << " " << meshDim;

    //meshLoader(mainElement, frontierElement, gaussType, meshDim); // Initialization of all quantities required.

    gmsh::model::list(modelNames);

    mainView.name = "MainView";
    mainView.tag = gmsh::view::add(mainView.name);
    mainView.dataType = "ElementNodeData";
    
    mainView.modelName = modelNames[0];
    mainView.data.resize(mainElement.elementTag.size(), std::vector<double>(mainElement.numNodes));

    solver(mainElement, frontierElement, mainView, simTime, timeInc, solvType, registration, debug); // Solving of the PDE with DG-FEM.

    gmsh::finalize(); // Closes gmsh
    return 0;
}
