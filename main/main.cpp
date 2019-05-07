#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "meshing.hpp"
#include "parameters.hpp"
#include "solver.hpp"
#include "structures.hpp"

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

    if(simulation.uNum == 6)
    {
        View EView; // View of the results.
        View HView; // View of the results.

        EView.name = "EView";
        EView.tag = gmsh::view::add(EView.name);
        EView.dataType = "ElementNodeData";

        HView.name = "HView";
        HView.tag = gmsh::view::add(HView.name);
        HView.dataType = "ElementNodeData";
        
        EView.modelName = HView.modelName = modelNames[0];

        EView.data.resize(mainElement.elementTag.size(), std::vector<double>(3 * mainElement.numNodes));
        HView.data.resize(mainElement.elementTag.size(), std::vector<double>(3 * mainElement.numNodes));

        solver(mainElement, frontierElement, physicalGroups, EView, HView, simulation); // Solving of the PDE with DG-FEM.
    }

    else if(simulation.uNum == 1)
    {
        View view, bin; // View of the results.

        view.name = "View";
        view.tag = gmsh::view::add(view.name);
        view.dataType = "ElementNodeData";
        
        view.modelName = modelNames[0];

        view.data.resize(mainElement.elementTag.size(), std::vector<double>(mainElement.numNodes));

        solver(mainElement, frontierElement, physicalGroups, view, bin, simulation); // Solving of the PDE with DG-FEM.
    }

    gmsh::finalize(); // Closes gmsh
    return 0;
}
