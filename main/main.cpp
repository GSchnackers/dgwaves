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

    Element mainElement; // The main elements of the mesh.
    Element frontierElement; // The frontier elements of the mesh.

    View mainView; // View of the results.

    gmsh::initialize(argc, argv); // Initialization of gmsh library.
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::open(argv[1]);                            // reads the msh file

    mainView.name = "Main View";
    mainView.tag = gmsh::view::add(mainView.name);
    mainView.dataType = "ElementNodeData";

    meshLoader(mainElement, frontierElement); // Initialization of all quantities required.

    /* double t; // t is the time.
    std::size_t i, j;

    std::vector<double> u(mainElement.nodeTags.size(), 0); // vector of the quantity of interest u.
    std::vector<double> uGp(mainElement.elementTag.size() * mainElement.numGp, 0); // vector of u at the Gauss Points.
    std::vector<double> uNext(mainElement.nodeTags.size()); // Next nodal value.

    std::vector<std::vector<double>> data(u.size(), std::vector<double>(1)); // puts u in a more suitable form for addition to data analysis.

    for(i = 0; i < u.size(); ++i) data[i][0] = u[i]; // Passing the u vector in the data vector.
    

    gmsh::view::addModelData( mainView.tag, 0, mainView.name, mainView.dataType, mainElement.nodeTags, data, 0, 1);

    std::vector<double> fluxNumGp(mainElement.nodeTags.size() * mainElement.numGp * 3, 0); // Vector of numerical fluxes at the Gauss points in parametric coordinates.
    std::vector<double> fluxPhysGp(fluxNumGp.size()); // Vector of physical fluxes at the Gauss Points in parametric coordinates.
    std::vector<double> fluxPhys(mainElement.nodeTags.size() * 3); // Vector of physical fluxes at the nodes in real coordinates.

    std::vector<double> velocity(3 * mainElement.nodeTags.size(), 0);
    std::vector<double> velocityGp(3 * mainElement.elementTag.size() * mainElement.numGp, 0);
    std::vector<double> velocityGpFrontier(3 * mainElement.elementTag.size() * mainElement.numGp, 0);

    // Velocity initializer.
    for(i = 0; i < velocity.size(); ++i)
        if(!(i % 3)) velocity[i] = 1;

    for(i = 0; i < velocityGp.size(); ++i)
        if(!(i % 3)) velocityGp[i] = 1;

    for(i = 0; i < velocityGpFrontier.size(); ++i)
        if(!(i % 3)) velocityGpFrontier[i] = 1;

    for(t = 0; t = 1; t += 0.01){

        valGp(u, mainElement, uGp);
        fluxComp(mainElement, u, uGp, velocity, velocityGp, fluxPhysGp, fluxPhys);

    } */


    gmsh::finalize(); // Closes gmsh
    return 0;
}
