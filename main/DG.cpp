#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(argv[1]);

    std::vector<int> eleTypes;
    gmsh::model::mesh::getElementTypes(eleTypes, 2);
    if (eleTypes.size() != 1)
    {
        gmsh::logger::write("Hybrid meshes not handled in this example!",
                            "error");
        return 1;
    }

    // 2D elements
    std::vector<std::pair<int, int>> entities2D;
    gmsh::model::getEntities(entities2D, 2);
    
    int eleType2D = eleTypes[0];
    std::string name;
    int dim, order, numNodes2D;
    std::vector<double> paramCoord;
    gmsh::model::mesh::getElementProperties(eleType2D, name, dim, order,
                                            numNodes2D, paramCoord);

    int s2D = entities2D[0].second;
    std::vector<int> elementTags2D, nodeTags2D;
    gmsh::model::mesh::getElementsByType(eleType2D, elementTags2D, nodeTags2D, s2D);

    std::vector<double> intpts2D, bf2D;
    int numComp2D;
    gmsh::model::mesh::getBasisFunctions(eleType2D, "Gauss4", "IsoParametric",
                                         intpts2D, numComp2D, bf2D);

    std::vector<double> jac2D, det2D, pts2D;
    gmsh::model::mesh::getJacobians(eleType2D, "Gauss4", jac2D, det2D, pts2D, s2D);

    std::vector<double> functionM;
    int numElements2D = elementTags2D.size();
    int numGaussPoints2D = intpts2D.size()/4;
    
    for(std::size_t e = 0; e < numElements2D; e++)
        for(std::size_t i = 0; i < numNodes2D; i++)
            for(std::size_t j = 0; j < numNodes2D; j++)
            {
                for(std::size_t g = 0; g < numGaussPoints2D; g++)
                    functionM.push_back(bf2D[numNodes2D*g + i] * bf2D[numNodes2D*g + j]);
            }
    
    std::vector<double> matrixM;
    gaussIntegration(intpts2D, functionM, det2D, matrixM, numElements2D, numGaussPoints2D, numNodes2D);



    // 1D elements
    for (std::size_t i = 0; i < entities2D.size(); i++)
    {
        s2D = entities2D[i].second;

        // get the nodes on the edges of the 2D elements
        std::vector<int> nodes;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodes, s2D);

        // Sorting duplicates
        sorting(nodes);

        // Computation of the normals to the elements.
        std::vector<double> normal2D(nodes.size());
        normal(nodes, normal2D);

        int c = gmsh::model::addDiscreteEntity(1);

        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodes);
        neighbours(nodeTags2D, numNodes2D, elementTags2D, nodes);
    }
    
    gmsh::model::mesh::getElementTypes(eleTypes, 1);

    int eleType1D = eleTypes[0];
    std::vector<double> intpts1D, bf1D;
    int numComp1D;
    gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss3", "IsoParametric",
                                         intpts1D, numComp1D, bf1D);

    std::vector<std::pair<int, int>> entities1D;
    gmsh::model::getEntities(entities1D, 1);

    // entities1D.size() = 5 -> 4 lines and 1 lines for surface
    for (std::size_t i = 0; i < entities1D.size(); i++)
    {
        int c = entities1D[i].second;
        std::vector<int> elementTags1D, nodeTags1D;
        gmsh::model::mesh::getElementsByType(eleType1D, elementTags1D, nodeTags1D, c);
        std::vector<double> jac1D, det1D, pts1D;
        gmsh::model::mesh::getJacobians(eleType1D, "Gauss3", jac1D, det1D, pts1D, c);
    }

    gmsh::finalize();
    return 0;
}