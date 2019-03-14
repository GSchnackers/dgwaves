#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

int main(int argc, char **argv)
{
    // Check of arguments
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    // Initialization of Gmsh
    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(argv[1]);

    // Get types of 2D elements and check if the mesh is not hybrid
    std::vector<int> eleTypes;
    gmsh::model::mesh::getElementTypes(eleTypes, 2);
    if (eleTypes.size() != 1)
    {
        gmsh::logger::write("Hybrid meshes not handled in this example!",
                            "error");
        return 1;
    }

    // 2D elements
    // Get the different surfaces of the mesh
    std::vector<std::pair<int, int>> entities2D;
    gmsh::model::getEntities(entities2D, 2);

    // Properties of the 2D elements in the mesh
    int eleType2D = eleTypes[0];
    std::string name;
    int dim, order, numNodes2D;
    std::vector<double> paramCoord;
    gmsh::model::mesh::getElementProperties(eleType2D, name, dim, order,
                                            numNodes2D, paramCoord);

    // tag of entity of edges
    int c = 0;

    // tag and nodes of edges
    std::vector<int> tagElement1D;
    std::vector<int> edgeNodes1D;

    // Loop on surfaces in the mesh
    for (std::size_t i = 0; i < entities2D.size(); i++)
    {
        int s2D = entities2D[i].second;

        // Get 2D elements of type eleType2D
        std::vector<int> elementTags2D, nodeTags2D;
        gmsh::model::mesh::getElementsByType(eleType2D, elementTags2D, nodeTags2D, s2D);

        // Get basis functions of 2D elements
        std::vector<double> intpts2D, bf2D;
        int numComp2D;
        gmsh::model::mesh::getBasisFunctions(eleType2D, "Gauss4", "IsoParametric",
                                            intpts2D, numComp2D, bf2D);

        // Get jacobian and its determinant of 2D elements
        std::vector<double> jac2D, det2D, pts2D;
        gmsh::model::mesh::getJacobians(eleType2D, "Gauss4", jac2D, det2D, pts2D, s2D);




        // function to integrate with Gauss integration to get the matrix M
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get the nodes on the edges of the 2D elements
        std::vector<int> edgeNodes2D;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, edgeNodes2D, s2D);

        //list of nodes for each element : nodeTags2D

        //list of nodal values
        std::vector<double> u(nodeTags2D.size());

        //déclaration coordonnées
        std::vector<double> nodeCoord(3);
        std::vector<double> nodeCoordParam(3);

        double value(nodeTags2D.size());
        //initial condition
        for(std::size_t i=0; i<nodeTags2D.size(); i++){
            gmsh::model::mesh::getNode(nodeTags2D[i], nodeCoord, nodeCoordParam);
            initialCondition(nodeCoord, value);
            u[i]=value;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////

        

        // Add an entity to contain the sorted edges
        c = gmsh::model::addDiscreteEntity(1);
        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, tagElement1D, edgeNodes1D);


    }
    
    // Get type of 1D elements 
    gmsh::model::mesh::getElementTypes(eleTypes, 1);

    // Get basis functions of 1D elements
    int eleType1D = eleTypes[0];
    std::vector<double> intpts1D, bf1D;
    int numComp1D;
    gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss3", "IsoParametric",
                                         intpts1D, numComp1D, bf1D);

    // Get jacobian and its determinant of 1D elements
    std::vector<double> jac1D, det1D, pts1D;
    gmsh::model::mesh::getJacobians(eleType1D, "Gauss3", jac1D, det1D, pts1D, c);


    ///////////////////////////////////////////////////////////////////////////////////////////////////


     int NumNodesSide = edgeNodes1D.size()/tagElement1D.size(); // number of nodes per side 
     int NumGaussPoint1D = det1D.size()/tagElement1D.size(); // number of gauss point per side


     /////////////////////////////////////////////////////////////////////////////////
     //trier edgeNodes1D , tagElement1D , det1D pour ne plus qu'il y ai de doublon
     //les vecteurs trié seront edgeNodes1DSorted , tagElement1DSorted , det1DSorted
     /////////////////////////////////////////////////////////////////////////////////

    std::vector<int> edgeNodes1DSorted;
    std::vector<int> tagElement1DSorted;
    std::vector<double> det1DSorted;

    //nodes of the first element
    for(size_t i=0; i<NumNodesSide; i++){

        edgeNodes1DSorted.push_back(edgeNodes1D[i]);
     }

    //tag of the first element
     tagElement1DSorted.push_back(tagElement1D[0]);

    //Jacobians of the first element
    for(size_t i=0; i<NumGaussPoint1D; i++){

        det1DSorted.push_back(det1D[i]);
     }
    
    //les indices i et j sont les indices du premier noeud de chaque edge    
     for(size_t i=0; i<edgeNodes1D.size(); i+= NumNodesSide){

         for(size_t j=0; j<edgeNodes1DSorted.size(); j+= NumNodesSide){

             // Check if edges is already in sortingNodes in the same direction
            if(edgeNodes1D[i] == edgeNodes1DSorted[j] && edgeNodes1D[i+ NumNodesSide -1] == edgeNodes1DSorted[j+ NumNodesSide -1])
            {
                break;
            }
            // Check if edges is already in sortingNodes in the opposite direction
            else if(edgeNodes1D[i] == edgeNodes1DSorted[j+ NumNodesSide -1] && edgeNodes1D[i+ NumNodesSide -1] == edgeNodes1DSorted[j])
            {
                break;
            }

            // If the edge is not already in sortingNodes, we add it
            if(j+NumNodesSide == edgeNodes1DSorted.size())
            {
                //fill edgeNodes1DSorted
                for(size_t n=0; n<NumNodesSide; n++){

                    edgeNodes1DSorted.push_back(edgeNodes1D[i+n]);
                }

                //fill tagElement1DSorted
                tagElement1DSorted.push_back(tagElement1D[i/NumNodesSide]);

                //fill det1DSorted
                for(size_t g=0; g<NumGaussPoint1D; g++){

                    det1DSorted.push_back(det1D[(i/NumNodesSide)*NumGaussPoint1D + g]);
                    //(i/NumNodesSide) est le numéro de l'edge 
                    //(i/NumNodesSide)*NumGausspoint1D + g est l'indice du point de gauss de l'edge

                }

            }//fin du if()

         }// fin de boucle sur j   
     }//fin de boucle sur i


/////////////////////////////////////////////
//Calculer les normales des éléments triés//
////////////////////////////////////////////
std::vector<double> normale (tagElement1DSorted*2); // car 2 composantes par edge

std::vector<double> nodeCoord1, nodeCoordParam1;
std::vector<double> nodeCoord2, nodeCoordParam2;

for(size_t i=0; i<edgeNodes1DSorted.size(); i+= NumNodesSide){

//calcul de la normale au bord
gmsh::model::mesh::getNode(edgeNodes1DSorted[i], nodeCoord1, nodeCoordParam1);
gmsh::model::mesh::getNode(edgeNodes1DSorted[i+NumNodesSide-1], nodeCoord2, nodeCoordParam2);

// Computation of the normal. n = (-y , x)/(x^2+y^2)^(1/2)

normale[(i/NumNodesSide)*2] = nodeCoord1[1] - nodeCoord2[1]; // -y
normale[(i/NumNodesSide)*2+1] = nodeCoord2[0] - nodeCoord1[0]; // x
//(i/NumNodesSide)*2 est la place dans normale seulement le edge

double norm = sqrt(normale[(i/NumNodesSide)*2] * normale[(i/NumNodesSide)*2] + \
                     normale[(i/NumNodesSide)*2+1] * normale[(i/NumNodesSide)*2+1]); //sqrt(x^2 + y^2)

// Final norm.
normale[(i/NumNodesSide)*2] /= norm;
n_2 /= norm;

}










/*
        // Get the neighbourhood of 2D elements
        std::vector<int> neighbourhood(nodes.size());
        neighbours(nodeTags2D, numNodes2D, elementTags2D, nodes, neighbourhood);

        // Computation of the normals to the elements.
        std::vector<double> normal2D(nodes.size());
        normal(nodes, normal2D);

*/

    gmsh::finalize();
    return 0;
}