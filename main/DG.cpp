#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "functions.h"

int main(int argc, char **argv)
{
    // coefF = [a_x ; a_y] the speed of the transport
    std::vector<double> coefF(2);

    // The user has to choose the values he wants for coefF
    coefF[0] = 1; //example
    coefF[1] = 0; //example

///////////////////////////////////////////////////////////////////////////////////////////////////////

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

    std::vector<std::string> names;
    gmsh::model::list(names);

    // Create a new post-processing view
    int viewtag = gmsh::view::add("my results");

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
    
    // Get 2D elements of type eleType2D
    std::vector<int> elementTags2D, nodeTags2D;

    // tag and nodes of edges
    std::vector<int> tagElement1D;
    std::vector<int> edgeNodes1D;

    // get the nodes on the edges of the 2D elements
    std::vector<int> edgeNodes2D;
    

    // Loop on surfaces in the mesh
    //for (std::size_t i = 0; i < entities2D.size(); i++)
    //{
        int s2D = entities2D[0].second;

        // Get 2D elements of type eleType2D
        gmsh::model::mesh::getElementsByType(eleType2D, elementTags2D, nodeTags2D, s2D);

        // Get basis functions of 2D elements
        std::vector<double> intpts2D, bf2D, gradIntPts2D, gradbf2D;
        int numComp2D, gradNumComp2D;
        gmsh::model::mesh::getBasisFunctions(eleType2D, "Gauss4", "IsoParametric",
                                            intpts2D, numComp2D, bf2D);

        // getBasisFunctions with option "GradLagrange" puts 
        // [g1 df1/du, g1 df1/dv, g1 df1/dw, ..., g1 dfC/dw, g2 df1/du, ...] in basis function
        gmsh::model::mesh::getBasisFunctions(eleType2D, "Gauss4", "GradLagrange",
                                            gradIntPts2D, gradNumComp2D, gradbf2D);
        
        // Get jacobian and its determinant of 2D elements
        std::vector<double> jac2D, det2D, pts2D, jac2DInverted;
        gmsh::model::mesh::getJacobians(eleType2D, "Gauss4", jac2D, det2D, pts2D, s2D);

        std::vector<double> detS2D(det2D.size());
        

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
                    for(std::size_t g = 0; g < numGaussPoints2D; g++){
                        functionM.push_back(bf2D[numNodes2D*g + i] * bf2D[numNodes2D*g + j]);
                    }
        
        std::vector<double> matrixM, matrixM_Inverted;
        gaussIntegration(intpts2D, functionM, det2D, matrixM, numElements2D, numGaussPoints2D, numNodes2D);

        std::vector<double> matrix_tmp(numNodes2D*numNodes2D), matrix_tmpInverted(numNodes2D*numNodes2D);

        // loop on 2D element and putting the matrices which are concatenated in a matrix for one element : matrix_tmp
        for(std::size_t e=0; e<elementTags2D.size(); e++){

            for(std::size_t i=0; i<numNodes2D; i++){
                for(std::size_t j=0; j<numNodes2D; j++){
                    matrix_tmp[i*numNodes2D + j] = matrixM[e*numNodes2D*numNodes2D + i*numNodes2D + j];
                }
            }
            invert(matrix_tmp, matrix_tmpInverted);

            matrixM_Inverted.insert( matrixM_Inverted.end(), matrix_tmpInverted.begin(), matrix_tmpInverted.end() );
        
        }



////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////     Matrix S     //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // attention la matrice jacobienne est toujours 3D
        std::vector<double> matrix3x3_tmp(9), matrix3x3_tmpInverted(9);
        // loop on 2D element and putting the matrices which are concatenated in a matrix for one element : matrix3x3_tmp
        for(std::size_t e=0; e<elementTags2D.size();e++){
            for(std::size_t g=0; g<numGaussPoints2D; g++){

                for(std::size_t i=0; i<3; i++){
                    for(std::size_t j=0; j<3; j++){
                        matrix3x3_tmp[i*3 + j] = jac2D[e*numGaussPoints2D*9 + g*9 + i*3 + j];
                    }
                }
                invert(matrix3x3_tmp, matrix3x3_tmpInverted);
                
                jac2DInverted.insert(jac2DInverted.end(), matrix3x3_tmpInverted.begin(), matrix3x3_tmpInverted.end());
                detS2D[numGaussPoints2D*e + g] = 1;
            }
        }// fin d'invertion
        /*
        for(std::size_t e=0; e<elementTags2D.size();e++){
            for(std::size_t g=0; g<numGaussPoints2D; g++){
                std::cout << "e" << std::to_string(e) << " g" << std::to_string(g)\
                << " x = " << std::to_string(intpts2D[4*g]) << " | y = " << std::to_string(intpts2D[4*g+1]) << "\n";
                for(std::size_t i=0; i<3; i++){
                    for(std::size_t j=0; j<3; j++){
                        std::cout << std::to_string(jac2D[e*numGaussPoints2D*9 + g*9 + i*3 + j]);
                    }
                    std::cout << "\n";
                }
                std::cout << "\n";
            }
        }
        */

        // function to integrate with Gauss integration to get the matrix S
        std::vector<double> functionS;
        double dfdx = 0, dfdy = 0;

        // [g1 df1/du, g1 df1/dv, g1 df1/dw, g1 df2/du ..., g1 dfN/dw, g2 df1/du, ...] in gradbf2D (N number of nodes of 2D element)
        
        for(std::size_t e = 0; e < numElements2D; e++){
            for(std::size_t i = 0; i < numNodes2D; i++)
                for(std::size_t j = 0; j < numNodes2D; j++)
                    for(std::size_t g = 0; g < numGaussPoints2D; g++){
                        for(std::size_t k = 0; k < 2; k++){
                            dfdx += gradbf2D[3*numNodes2D*g + 3*j + k] * jac2DInverted[e*9*numGaussPoints2D + 9*g + k];
                            dfdy += gradbf2D[3*numNodes2D*g + 3*j + k] * jac2DInverted[e*9*numGaussPoints2D + 9*g + k + 3];
                        }

                        functionS.push_back((coefF[0]*dfdx + coefF[1]*dfdy)*bf2D[numNodes2D*g + i]);
                        dfdx = 0;
                        dfdy = 0;
                    }
        }
        
        std::vector<double> matrixS;
        gaussIntegration(intpts2D, functionS, det2D, matrixS, numElements2D, numGaussPoints2D, numNodes2D);
        
        for(std::size_t e = 0; e < numElements2D; e++){
            for(std::size_t i = 0; i < numNodes2D; i++){
                for(std::size_t j = 0; j < numNodes2D; j++)
                {
                    std::cout << std::to_string(matrixS[numNodes2D*(numNodes2D*e + i) + j]) << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// Nodal values of u and list of nodes for each element ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // get the nodes on the edges of the 2D elements
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, edgeNodes2D, s2D);

        //list of nodes for each element : nodeTags2D

        //list of nodal values
        std::vector<double> u(nodeTags2D.size());

        //list of du/dt
        std::vector<double> dudt(u.size());

        //all data
        std::vector<double> tmp(numNodes2D);
        std::vector<std::vector<double>> data;

        //déclaration coordonnées
        std::vector<double> nodeCoord(3);
        std::vector<double> nodeCoordParam(3);

        double value;
        //initial condition
        for(std::size_t i=0; i<nodeTags2D.size(); i++){
            gmsh::model::mesh::getNode(nodeTags2D[i], nodeCoord, nodeCoordParam);
            initialCondition(nodeCoord, value);
            u[i]=value;
        }

        for(size_t i = 0; i < u.size(); i++){
            std::cout << "u[" << std::to_string(i) << "] : " << std::to_string(u[i]) << "\n";
        }

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Entity 1D ////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
        
        // Add an entity to contain the sorted edges
        c = gmsh::model::addDiscreteEntity(1);
        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, edgeNodes2D);

    //}
    
    // Get type of 1D elements 
    gmsh::model::mesh::getElementTypes(eleTypes, 1);

    // Properties of the 2D elements in the mesh
    eleType1D = eleTypes[0];
    std::string name1D;
    int dim1D, order1D, numNodes1D;
    std::vector<double> paramCoord1D;
    gmsh::model::mesh::getElementProperties(eleType1D, name1D, dim1D, order1D,
                                            numNodes1D, paramCoord1D);

    // Get basis functions of 1D elements
    std::vector<double> intpts1D, bf1D;
    int numComp1D;
    gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss3", "IsoParametric",
                                         intpts1D, numComp1D, bf1D);

    // Get jacobian and its determinant of 1D elements
    std::vector<double> jac1D, det1D, pts1D;
    gmsh::model::mesh::getJacobians(eleType1D, "Gauss3", jac1D, det1D, pts1D, c);

    // Get 1D elements of type eleType1D
    std::vector<int> elementTags1D, nodeTags1D;
    gmsh::model::mesh::getElementsByType(eleType1D, elementTags1D, nodeTags1D, c);



    ///////////////////////////////////////////////////////////////////////////////////////////////////

    tagElement1D = elementTags1D;
    edgeNodes1D = nodeTags1D;

    int NumNodesSide = edgeNodes1D.size()/tagElement1D.size(); // number of nodes per side 
    int NumGaussPoint1D = det1D.size()/tagElement1D.size(); // number of gauss point per side

    std::cout << "NumNodesSide = " << std::to_string(NumNodesSide) << "\n";
    std::cout << "NumGaussPoint1D = " << std::to_string(NumGaussPoint1D) << "\n";


    ///////////////////////////////////////////////////////////////////////////////////////
    //// trier, pour ne plus qu'il y ai de doublon :
    //// edgeNodes1D , tagElement1D , det1D    
    ////
    //// les vecteurs trié seront :
    //// edgeNodes1DSorted , tagElement1DSorted , det1DSorted
    ///////////////////////////////////////////////////////////////////////////////////////

    std::vector<int> edgeNodes1DSorted;
    std::vector<int> tagElement1DSorted;
    std::vector<double> det1DSorted;

    //nodes of the first element
    for(std::size_t i=0; i<NumNodesSide; i++){

        edgeNodes1DSorted.push_back(edgeNodes1D[i]);
    }

    //tag of the first element
     tagElement1DSorted.push_back(tagElement1D[0]);

    //Jacobians of the first element
    for(std::size_t i=0; i<NumGaussPoint1D; i++){

        det1DSorted.push_back(det1D[i]);
    }
    
    //les indices i et j sont les indices du premier noeud de chaque edge    
    for(std::size_t i=0; i<edgeNodes1D.size(); i+= NumNodesSide){

        for(std::size_t j=0; j<edgeNodes1DSorted.size(); j+= NumNodesSide){

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
                for(std::size_t n=0; n<NumNodesSide; n++){

                    edgeNodes1DSorted.push_back(edgeNodes1D[i+n]);
                }

                //fill tagElement1DSorted
                tagElement1DSorted.push_back(tagElement1D[i/NumNodesSide]);

                //fill det1DSorted
                for(std::size_t g=0; g<NumGaussPoint1D; g++){

                    det1DSorted.push_back(det1D[(i/NumNodesSide)*NumGaussPoint1D + g]);
                    //(i/NumNodesSide) est le numéro de l'edge 
                    //(i/NumNodesSide)*NumGausspoint1D + g est l'indice du point de gauss de l'edge

                }

            }//fin du if()

        }// fin de boucle sur j   
    }//fin de boucle sur i

    /* TEST
    std::cout << "edgeNodes1DSorted.size() = " << std::to_string(edgeNodes1DSorted.size()) << "\n";
    std::cout << "tagElement1DSorted.size() = " << std::to_string(tagElement1DSorted.size()) << "\n";
    std::cout << "det1DSorted.size() = " << std::to_string(det1DSorted.size()) << "\n";
    */

   // TEST
    for(size_t i = 0; i < edgeNodes1DSorted.size(); i++){
        std::cout << "edgeNodes1DSorted[" << std::to_string(i) << "]  : " << std::to_string(edgeNodes1DSorted[i]) << "\n";
    }

    //////////////////////////////////////////////////////////////////////////////
    ///////////////// Calculer les normales des éléments triés ///////////////////
    //////////////////////////////////////////////////////////////////////////////

    std::vector<double> normal(tagElement1DSorted.size()*2); // car 2 composantes par edge

    std::vector<double> nodeCoord1, nodeCoordParam1;
    std::vector<double> nodeCoord2, nodeCoordParam2;

    for(std::size_t i=0; i<edgeNodes1DSorted.size(); i+= NumNodesSide){

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

    //Test
    /*
    for(size_t i = 0; i < normal.size(); i++){
        std::cout << "normal[" << std::to_string(i) << "]  : " << std::to_string(normal[i]) << "\n";
    }
    */

    ////////////////////////////////////////////////////////////////////////////
    ////////// Find the neighbours of the edges and Find the BC tags ///////////
    ////////////////////////////////////////////////////////////////////////////
    int neighbour1D_tmp;
    std::vector<int> neighbours1D(tagElement1DSorted.size()*2,-1); // pour la taille "2*"car 2 voisins par edge,
                                                                // on initialise à -1 ainsi nous sauront lorsqu'un côté n'a pas de voisin

    //nombre de côtés de l'element 2D = nombre de noeuds sur les côtés divisé par (nombre de noeud par côtés)
    int NumSide2D = edgeNodes2D.size()/(NumNodesSide*elementTags2D.size());

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
    for(std::size_t i=0; i<edgeNodes1DSorted.size(); i+= NumNodesSide){

        for(std::size_t j=0; j<edgeNodes2D.size(); j+= NumNodesSide){
            // Check if the edge is common to the edge of one 2D element
            if((edgeNodes2D[j] == edgeNodes1DSorted[i] && edgeNodes2D[j+ NumNodesSide -1] == edgeNodes1DSorted[i+ NumNodesSide -1]) \
                ||(edgeNodes2D[j] == edgeNodes1DSorted[i+ NumNodesSide -1] && edgeNodes2D[j+ NumNodesSide -1] == edgeNodes1DSorted[i])){
                
                neighbour1D_tmp = j/(NumSide2D*NumNodesSide);
                // (j/NumNodesSide) is the number("index") of the 2D element

                // look at the inner product between the vector constructed with two node coordinates of the 2D element (one positioned on the edge)
                // and the vector normal to the edge to know if the normal is in the conventional direction or not

                // index_tmp contient l'index du noeud suivant de l'élément 2D 
                index_tmp = (j+ 2*NumNodesSide) % (NumSide2D*NumNodesSide);

                //get the coordinates of the begin and end nodes of the next edge of the 2D element
                gmsh::model::mesh::getNode(edgeNodes2D[j+ NumNodesSide -1], nodeCoord1, nodeCoordParam1);
                gmsh::model::mesh::getNode(edgeNodes2D[neighbour1D_tmp*NumNodesSide*NumSide2D + index_tmp], nodeCoord2, nodeCoordParam2);

                //inner product
                innerProduct = normal[i/(NumNodesSide)*2]*(nodeCoord2[0]-nodeCoord1[0]) + \
                                normal[i/(NumNodesSide)*2+1]*(nodeCoord2[1]-nodeCoord1[1]);

                //if the neighbour is positionned conventionnaly with respect to the normal, its number is registered in first position,
                //if not, its number is registered in second position.
                if(innerProduct <= 0){
                    neighbours1D[2*i/(NumNodesSide)] = neighbour1D_tmp;
                }
                else{
                    neighbours1D[2*i/(NumNodesSide) + 1] = neighbour1D_tmp;
                }

            } // end check if neighbour

        }// fin de boucle sur j   
    

        //fill the tags for BC in nodeTags2DPlusBC
        if(neighbours1D[i/(NumNodesSide)*2] == -1 || neighbours1D[i/(NumNodesSide)*2 + 1] == -1){
            
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

    for(size_t i = 0; i < neighbours1D.size(); i++){
        std::cout << "neighbours1D[" << std::to_string(i) << "]  : " << std::to_string(neighbours1D[i]) << "\n";
    }

    // TEST
    /*
    for(size_t i = 0; i < edgeNodes2D.size(); i++){
        std::cout << "edgeNodes2D[" << std::to_string(i) << "]  : " << std::to_string(edgeNodes2D[i]) << "\n";
    }
    */
    for(size_t i = 0; i < nodeTags2DPlusBC.size(); i++){
        std::cout << "nodeTags2DPlusBC[" << std::to_string(i) << "]  : " << std::to_string(nodeTags2DPlusBC[i]) << "\n";
    }
    

    ////////////////////////////////////////////////////////////////////////////////////
    // association of the nodes of edgeNodes1DSorted with their indices in nodeTags2D //
    ////////////////////////////////////////////////////////////////////////////////////

    std::vector<int> indicesNei1(edgeNodes1DSorted.size(),-1);
    std::vector<int> indicesNei2(edgeNodes1DSorted.size(),-1);

    for(size_t i=0; i<edgeNodes1DSorted.size(); i++){

        // First neighbour
        // check if there is a neighbour
        if(neighbours1D[2*(i/NumNodesSide)] != -1){
            // as we know the number of the neighbour element, we will only loop on its nodes in the general list : nodeTags2D
            for(size_t j=neighbours1D[2*(i/NumNodesSide)]*numNodes2D; j<(neighbours1D[2*(i/NumNodesSide)] + 1)*numNodes2D; j++){

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
        if(neighbours1D[2*(i/NumNodesSide) + 1] != -1){
            for(size_t j=neighbours1D[2*(i/NumNodesSide) + 1]*numNodes2D; j<(neighbours1D[2*(i/NumNodesSide) + 1] + 1)*numNodes2D; j++){

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

    for(size_t i = 0; i < indicesNei1.size(); i++){
        std::cout << "indicesNei1[" << std::to_string(i) << "]  : " << std::to_string(indicesNei1[i]) << "\n";
    }

    for(size_t i = 0; i < indicesNei2.size(); i++){
        std::cout << "indicesNei2[" << std::to_string(i) << "]  : " << std::to_string(indicesNei2[i]) << "\n";
    }


    ///////////////////////////////////////////////////////////////////////
    //////////////////////// upwind variable //////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    /* 
    upwind register if the inner product of the normal and coefF (the speed of transport) is positive or negative
    This information will be used to know which neighbour is upwind
    - if upwind = 1 --> the first neighbour is upwind
    - if upwind = -1 --> the second neighbour is upwind
    */

    std::vector<int> upwind(tagElement1DSorted.size());

    for(std::size_t i=0; i<upwind.size(); i++){

        if(normal[i*2]*coefF[0] + normal[i*2+1]*coefF[1]>=0){
            upwind[i] = 1;
        }
        else{
            upwind[i] = -1;
        }
    }

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////// Matrix F ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    std::vector<double> matrixF(tagElement1DSorted.size()*NumNodesSide*NumNodesSide);

    //gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss3", "IsoParametric",
    //                                     intpts1D, numComp1D, bf1D);

    //gmsh::model::mesh::getJacobians(eleType1D, "Gauss3", jac1D, det1D, pts1D, c);  det1DSorted

    //loop for each edge
    for(std::size_t ed=0; ed < tagElement1DSorted.size(); ed++){
        //loop for i of F_{ij}
        for(std::size_t i=0; i < NumNodesSide; i++){
            //loop for j of F_{ij}
            for(std::size_t j=0; j < NumNodesSide; j++){
                //loop for each gauss point
                for(std::size_t g=0; g < NumGaussPoint1D; g++){

                    matrixF[ed*NumNodesSide*NumNodesSide + i*NumNodesSide + j] += \
                     (normal[ed*2]*coefF[0] + normal[ed*2+1]*coefF[1]) * bf1D[g*NumGaussPoint1D + i] * bf1D[g*NumGaussPoint1D + j] \
                     * intpts1D[3 + 4*g] * det1DSorted[ed*NumGaussPoint1D + g];

                }
            }
        }
    }

    for(std::size_t ed=0; ed < tagElement1DSorted.size(); ed++){
        std::cout << "e " << std::to_string(ed) << "\n";
        //loop for i of F_{ij}
        for(std::size_t i=0; i < NumNodesSide; i++){
            //loop for j of F_{ij}
            for(std::size_t j=0; j < NumNodesSide; j++){
                std::cout << std::to_string(matrixF[ed*NumNodesSide*NumNodesSide + i*NumNodesSide + j]) << " ";
            }
        std::cout << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "TIME LOOP\n";
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////           TIME LOOP           /////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double time = 0;
    double timeStep = 0.2;
    double endTime = 10;
    int numStep = endTime/timeStep;

    std::vector<double> Su(elementTags2D.size()*numNodes2D);

    std::string modelName = names[0];
    std::string dataType = "ElementNodeData";

    //fill data with u
    for(std::size_t e=0; e<elementTags2D.size(); e++){
        for(std::size_t i=0; i<numNodes2D; i++){
            tmp[i] = u[e*numNodes2D + i];
        }
        data.push_back(tmp);
    }

    gmsh::view::addModelData(viewtag, 0, modelName, dataType, elementTags2D, data, time, 1);

    // declaration vector F (time dependent)
    std::vector<double> vectorF(nodeTags2D.size());

    for(std::size_t step = 1; time < endTime; step++){

        // BC
        for(std::size_t i=nodeTags2D.size(); i<nodeTags2DPlusBC.size(); i++){
            gmsh::model::mesh::getNode(nodeTags2DPlusBC[i], nodeCoord, nodeCoordParam);
            boundaryConditions(nodeCoord, time, value);
            uPlusBC[i]=value;
        }

        /* TEST
        for(size_t i = 0; i < uPlusBC.size(); i++){
            std::cout << "uPlusBC[" << std::to_string(i) << "]  : " << std::to_string(uPlusBC[i]) << " at time " << std::to_string(time) << "\n";
        }
        */

        // initialisation to 0 of vector F
        for(std::size_t i=0; i<vectorF.size(); i++){
            vectorF[i] = 0;
        }

        // computation vector F
        for(std::size_t ed=0; ed<tagElement1DSorted.size(); ed++){
            if(upwind[ed] == 1){

                if(neighbours1D[ed*2] != -1){
                    //fill vectorF .... (produit mat)
                    //loop on the nodes of the edge 
                    for(std::size_t i=0; i<NumNodesSide; i++){
                        for(std::size_t j=0; j<NumNodesSide; j++){
                            vectorF[indicesNei1[ed*NumNodesSide + i]] += \
                            -(matrixF[ed*NumNodesSide*NumNodesSide + i*NumNodesSide + j] * uPlusBC[indicesNei1[ed*NumNodesSide + j]]);
                        }
                    }

                    if(neighbours1D[ed*2 + 1] != -1){
                        for(std::size_t copy=0; copy<NumNodesSide; copy++){
                            vectorF[indicesNei2[ed*NumNodesSide + copy]] += -vectorF[indicesNei1[ed*NumNodesSide + copy]];
                        }
                    }
                }
                if(neighbours1D[ed*2] == -1){
                    //vectorF[indicesNei2] = -... (produit mat)
                    for(std::size_t i=0; i<NumNodesSide; i++){
                        for(std::size_t j=0; j<NumNodesSide; j++){
                            vectorF[indicesNei2[ed*NumNodesSide + i]] += \
                            matrixF[ed*NumNodesSide*NumNodesSide + i*NumNodesSide + j] * uPlusBC[indicesNei1[ed*NumNodesSide + j]];
                        }
                    }

                }
            }
        
            if(upwind[ed] == -1){
            // same as upwind == 1 except we take "indicesNei2" to compute the flow            
                if(neighbours1D[ed*2] != -1){
                    //fill vectorF .... (produit mat)
                    //loop on the nodes of the edge 
                    for(std::size_t i=0; i<NumNodesSide; i++){
                        for(std::size_t j=0; j<NumNodesSide; j++){
                            vectorF[indicesNei1[ed*NumNodesSide + i]] += \
                            -(matrixF[ed*NumNodesSide*NumNodesSide + i*NumNodesSide + j] * uPlusBC[indicesNei2[ed*NumNodesSide + j]]);
                        }
                    }

                    if(neighbours1D[ed*2 + 1] != -1){
                        for(std::size_t copy=0; copy<NumNodesSide; copy++){
                            vectorF[indicesNei2[ed*NumNodesSide + copy]] += -vectorF[indicesNei1[ed*NumNodesSide + copy]];
                        }
                    }
                }
                if(neighbours1D[ed*2] == -1){
                    //vectorF[indicesNei2] = -... (produit mat)
                    for(std::size_t i=0; i<NumNodesSide; i++){
                        for(std::size_t j=0; j<NumNodesSide; j++){
                            vectorF[indicesNei2[ed*NumNodesSide + i]] += \
                            matrixF[ed*NumNodesSide*NumNodesSide + i*NumNodesSide + j] * uPlusBC[indicesNei2[ed*NumNodesSide + j]];
                        }
                    }

                }        
            }
        }// end computation vector F

        // Su à zéro
        for(std::size_t el=0; el<elementTags2D.size(); el++){
            for(std::size_t i=0; i<numNodes2D; i++){

                Su[el*numNodes2D + i] = 0;
            }
        }

        //Computation of Su = S.u
        for(std::size_t el=0; el<elementTags2D.size(); el++){
            for(std::size_t i=0; i<numNodes2D; i++){
                for(std::size_t j=0; j<numNodes2D; j++){

                    Su[el*numNodes2D + i] += matrixS[el*numNodes2D*numNodes2D + i*numNodes2D + j] * u[el*numNodes2D + j];
                }
            }
        }

        // dudt à zéro
        for(std::size_t el=0; el<elementTags2D.size(); el++){
            for(std::size_t i=0; i<numNodes2D; i++){

                dudt[el*numNodes2D + i] = 0;
            }
        }

        // du/dt = M^{-1} (S.u + F)
        for(std::size_t el=0; el<elementTags2D.size(); el++){
            for(std::size_t i=0; i<numNodes2D; i++){
                for(std::size_t j=0; j<numNodes2D; j++){

                    dudt[el*numNodes2D + i] += matrixM_Inverted[el*numNodes2D*numNodes2D + i*numNodes2D + j]* \
                                                    (Su[el*numNodes2D + j] + vectorF[el*numNodes2D + j]);
                }
            }
        }

        // Forward Euler method
        Forward_Euler_method(u, timeStep, dudt);

        //fill data with u
        for(std::size_t e=0; e<elementTags2D.size(); e++){

            for(std::size_t i=0; i<numNodes2D; i++){
                data[e][i] = u[e*numNodes2D + i];
            }
        }

        // Backup of u(t+dt)
        gmsh::view::addModelData(viewtag, step, modelName, dataType, elementTags2D, data, time, 1);

        time += timeStep;
    }

    gmsh::view::write(viewtag, std::string("results.msh"));


    gmsh::finalize();
    return 0;
}



/* Notes sur le code

ligne 370, peut être pas assez robuste, si jamais le produit scalaire vaut exactement 0,
 le deuxième voisin va être écrit en effaçant le premier:

 if(innerProduct >= 0){
                    neighbours1D[i/(NumNodesSide)] = neighbour1D_tmp;
                }
                else{
                    neighbours1D[i/(NumNodesSide) + 1] = neighbour1D_tmp;
                }



*/