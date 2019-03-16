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

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////     Matrix M     //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////




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
//////////////// Nodal values of u and list of nodes for each element ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // c = [a_x ; a_y] the speed of the transport
        std::vector<double> c(2);
        // The user has to choose the values he wants for c
        c[0] = 3.78; //example
        c[1] = 1.41; //example

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
////////////////////////////// Entity 1D ////////////////////////////////////////////////////
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


    ///////////////////////////////////////////////////////////////////////////////////////
    //// trier edgeNodes1D , tagElement1D , det1D pour ne plus qu'il y ai de doublon   ////
    //// les vecteurs trié seront edgeNodes1DSorted , tagElement1DSorted , det1DSorted ////
    ///////////////////////////////////////////////////////////////////////////////////////

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

            // If the edge is not already in edgeNodes1DSorted, we add it
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


    //////////////////////////////////////////////////////////////////////////////
    ///////////////// Calculer les normales des éléments triés ///////////////////
    //////////////////////////////////////////////////////////////////////////////

    std::vector<double> normal(tagElement1DSorted.size()*2); // car 2 composantes par edge

    std::vector<double> nodeCoord1, nodeCoordParam1;
    std::vector<double> nodeCoord2, nodeCoordParam2;

    for(size_t i=0; i<edgeNodes1DSorted.size(); i+= NumNodesSide){

    //calcul de la normale au bord
    gmsh::model::mesh::getNode(edgeNodes1DSorted[i], nodeCoord1, nodeCoordParam1);
    gmsh::model::mesh::getNode(edgeNodes1DSorted[i+NumNodesSide-1], nodeCoord2, nodeCoordParam2);

    // Computation of the normal. n = (-y , x)/(x^2+y^2)^(1/2)

    normal[(i/NumNodesSide)*2] = nodeCoord1[1] - nodeCoord2[1]; // -y
    normal[(i/NumNodesSide)*2+1] = nodeCoord2[0] - nodeCoord1[0]; // x
    //(i/NumNodesSide)*2 est la place dans normale seulement le edge

    double norm = sqrt(normal[(i/NumNodesSide)*2] * normal[(i/NumNodesSide)*2] + \
                        normal[(i/NumNodesSide)*2+1] * normal[(i/NumNodesSide)*2+1]); //sqrt(x^2 + y^2)

    // Final norm.
    normal[(i/NumNodesSide)*2] /= norm;
    normal[(i/NumNodesSide)*2+1] /= norm;

    }


    ////////////////////////////////////////////////////////////////////////////
    ////////// Find the neighbours of the edges and Find the BC tags ///////////
    ////////////////////////////////////////////////////////////////////////////

    int neighbour1D_tmp;
    std::vector<int> neighbours1D(tagElement1DSorted.size()*2,-1); // pour la taille "2*"car 2 voisins par edge,
                                                                // on initialise à -1 ainsi nous sauront lorsqu'un côté n'a pas de voisin

    //nombre de côtés de l'element 2D = nombre de noeuds sur les côtés divisé par (nombre de noeud par côtés)
    NumSide2D = edgeNodes2D.size()/(NumNodesSide*elementTags2D.size());

    int index_tmp;
    double innerProduct;

    // u concatenated with boundary conditions (we add BC later)
    std::vector<double> uPlusBC(u.size());

    // nodeTags2D concatenated with the tags of the nodes which have boundary conditions 
    //(we add those tags later, in the loop for the neighbours)
    std::vector<int> nodeTags2DPlusBC(nodeTags2D.size());

    for(std::size_t i=0; i<nodeTags2D.size(); i++){
        uPlusBC[i] = u[i];
        nodeTags2DPlusBC[i] = nodeTags2D[i];
    }

    int check1;
    int check2;


    //les indices i et j sont les indices du premier noeud de chaque edge    
    for(size_t i=0; i<edgeNodes1DSorted.size(); i+= NumNodesSide){

        for(size_t j=0; j<edgeNodes2D.size(); j+= NumNodesSide){

            // Check if the edge is common to the edge of one 2D element
            if((edgeNodes2D[j] == edgeNodes1DSorted[i] && edgeNodes2D[j+ NumNodesSide -1] == edgeNodes1DSorted[i+ NumNodesSide -1]) \
                ||(edgeNodes2D[j] == edgeNodes1DSorted[i+ NumNodesSide -1] && edgeNodes2D[j+ NumNodesSide -1] == edgeNodes1DSorted[i])){
                    
                neighbour1D_tmp = j/(NumSide2D*NumNodesSide);
                // (j/NumNodesSide) is the number("index") of the 2D element

                // look at the inner product between the vector constructed with two node coordinates of the 2D element (one positioned on the edge)
                // and the vector normal to the edge to know if the normal is in the conventional direction or not

                // index_tmp contient l'index du noeud suivant de l'élément 2D 
                index_tmp = j+ 2*NumNodesSide % (NumSide2D*NumNodesSide);

                //get the coordinates of the begin and end nodes of the next edge of the 2D element
                gmsh::model::mesh::getNode(edgeNodes2D[j+ NumNodesSide -1], nodeCoord1, nodeCoordParam1);
                gmsh::model::mesh::getNode(edgeNodes2D[neighbour1D_tmp*NumNodesSide*NumSide2D + index_tmp], nodeCoord2, nodeCoordParam2);

                //inner product
                innerProduct = normal[i/(NumNodesSide)]*(nodeCoord2[0]-nodeCoord1[0]) + \
                                normal[i/(NumNodesSide)+1]*(nodeCoord2[1]-nodeCoord1[1]);

                //if the neighbour is positionned conventionnaly with respect to the normal, its number is registered in first position,
                //if not, its number is registered in second position.
                if(innerProduct >= 0){
                    neighbours1D[i/(NumNodesSide)] = neighbour1D_tmp;
                }
                else{
                    neighbours1D[i/(NumNodesSide) + 1] = neighbour1D_tmp;
                }

            } // end check if neighbour

        }// fin de boucle sur j   

        //fill the tags for BC in nodeTags2DPlusBC
        if(neighbours1D[i/(NumNodesSide)] == -1 || neighbours1D[i/(NumNodesSide) + 1] == -1){
            
            check1 = 0;
            check2 = 0;
            // check if the tag of the node is already in the BC tags
            for(std::size_t c=nodeTags2D.size(); c<nodeTags2DPlusBC.size(); c++){
                if(nodeTags2DPlusBC[c] == edgeNodes1DSorted[i]){
                    check1 = 1;
                }
                if(nodeTags2DPlusBC[c] == edgeNodes1DSorted[i + NumNodesSide -1]){
                    check2 = 1;
                }
            }

            if(check1 == 0){
                nodeTags2DPlusBC.push_back(edgeNodes1DSorted[i]);
                uPlusBC.push_back(0);
            }
            if(check2 == 0){
                nodeTags2DPlusBC.push_back(edgeNodes1DSorted[i + NumNodesSide -1]);
                uPlusBC.push_back(0);
            }

        } // end fill tags of BC

    }//fin de boucle sur i

    ////////////////////////////////////////////////////////////////////////////////////
    // association of the nodes of edgeNodes1DSorted with their indices in nodeTags2D //
    ////////////////////////////////////////////////////////////////////////////////////

    std::vector<int> indicesNei1(edgeNodes1DSorted.size());
    std::vector<int> indicesNei2(edgeNodes1DSorted.size());

    for(size_t i=0; i<edgeNodes1DSorted.size(); i+= NumNodesSide){

        // First neighbour
        // check if there is a neighbour
        if(neighbours1D[i/(NumNodesSide)] != -1){
            // as we know the number of the neighbour element, we will only loop on its nodes in the general list : nodeTags2D
            for(size_t j=neighbours1D[i/(NumNodesSide)]*NumNodes2D; j<(neighbours1D[i/(NumNodesSide)] + 1)*NumNodes2D; j++){

                if(edgeNodes1DSorted[i] == nodeTags2D[j]){

                    indicesNei1[i] = j;
                }
            }
        }
        //if there is no neighbour
        else{
            // loop on the nodes of the Boundary Conditions
            for(size_t j=nodeTags2D.size(); j<nodeTags2DPlusBC.size(); j++){

                if(edgeNodes1DSorted[i] == nodeTags2DPlusBC[j]){

                    indicesNei1[i] = j;
                }
            }
        }

        // Second neighbour
        // Same steps as for the first neighbour
        if(neighbours1D[i/(NumNodesSide) + 1] != -1){
            for(size_t j=neighbours1D[i/(NumNodesSide) + 1]*NumNodes2D; j<(neighbours1D[i/(NumNodesSide) + 1] + 1)*NumNodes2D; j++){

                if(edgeNodes1DSorted[i] == nodeTags2D[j]){

                    indicesNei2[i] = j;
                }
            }
        }
        else{
            for(size_t j=nodeTags2D.size(); j<nodeTags2DPlusBC.size(); j++){

                if(edgeNodes1DSorted[i] == nodeTags2DPlusBC[j]){

                    indicesNei2[i] = j;
                }
            }
        }

    }// end loop on i


    ///////////////////////////////////////////////////////////////////////
    //////////////////////// upwind variable //////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    /* 
    upwind register if the inner product of the normal and c (the speed of transport) is positive or negative
    This information will be used to know which neighbour is upwind
    - if upwind = 1 --> the first neighbour is upwind
    - if upwind = -1 --> the second neighbour is upwind
    */

    std::vector<int> upwind(tagElement1D.size());

    for(std::size_t i=0; i<upwind.size(); i++){

        if(normal[i*2]*c[0] + normal[i*2+1]*c[1]>=0){
            upwind[i] = 1;
        }
        else{
            upwind[i] = -1;
        }
    }

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////// Matrix F ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////


























    gmsh::finalize();
    return 0;
}