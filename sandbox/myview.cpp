// This example generate some results over a mesh
//
// How to use this example?
//
// 1.build the code:
//     cd build
//     cmake -G "MinGW Makefiles" ..
//     cmake --build .      (or "mingw32-make")
// 2.generate a mesh from a geo file:
//     gmsh -2 -order 3 ..\sandbox\sea.geo
// 3.run the program with the msh as argument
//     bin\myview.exe ..\sandbox\sea.msh
// 4.display the generated msh data
//     gmsh data.msh





#include <gmsh.h>
#include <iostream>

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }
    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::open(argv[1]);                            // reads the msh file


    std::vector<std::string> names;
    gmsh::model::list(names);
    std::cout << "the file contains " << names.size() << " model names\n";


    // get all nodes
    std::vector<int> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoord);

    std::cout << "the mesh has " << nodeTags.size() << " nodeTags\n";
    std::cout << "the mesh has " << coords.size() << " coords\n";
    std::cout << "the mesh has " << parametricCoord.size() << " parametricCoord\n";

    // Create a new post-processing view
    int viewtag = gmsh::view::add("my results");

    std::vector<std::vector<double> > data(nodeTags.size());  // check this size!!

    int nstep = 20;
    double Lx = 20.;
    double v = 10.0;
    double tend = Lx/v;
    for(int step = 0; step < nstep; step++)
    {
        double time = tend*((double)step/nstep);

        for(int i=0; i<nodeTags.size(); ++i)
        {
            int tag = nodeTags[i];
            double x = coords[i*3+0];
            double y = coords[i*3+1];
            double z = coords[i*3+2];

            double val = sin(2*M_PI*(x-v*time)/Lx*2)+y;
            data[i].resize(1);
            data[i][0] = val;
        }

        std::string modelName = names[0];
        std::string dataType = "NodeData";

        gmsh::view::addModelData(viewtag, step, modelName, dataType,
                    nodeTags, data, time);
    }

    gmsh::view::write(viewtag, "data.msh");

    gmsh::finalize();
    return 0;
}
