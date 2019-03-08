// Implementation of the vector F

/* Input: 
- coef a_x and a_y
- values of u at the nodes (size = N_elements * N_basis_function)
- values of the basis function at gauss points of the edge (bfg = basis functions at gauss point) (size = N_elements * N_basis_function * N_gausspoint)
- nodes_unsorted (before mysorting - after getElementEdgeNodes - primary !!!!!!!!!!!!!!!!!!!!!!!!!!!!)
*/   

// Output: -  vector F (size N_elements * N_basis_function) attention, il faut l'initialiser à 0 après chaque pas de temps

#include <cstdio>
#include <iostream>
#include <vector>
#include "functions.h"
#include <gmsh.h>

void vectorF(std::vector<double> & vectorF, double a_x, double a_y, std::vector<double> u,
                  std::vector<double> bfg, std::vector<int> nodes_unsorted, int N_elements, int N_gauss_point){


    int N_basis_function = u.size()/N_elements;
    int N_nodes_per_element = nodes_unsorted.size()/(2*N_elements);
    int el_neighbour; //"index" of the neighboor element

    double nodeCoord1, nodeCoord2, nodeCoordParam1, nodeCoordParam2, nodeCoord3, nodeCoordParam3;
    int sign;


    //loop on all elements
    for(std::size_t el=0; el<N_element; el++){
        //loop on all edges of the element
        for(std::size_t n=0; n<N_nodes_per_element*2; n+=2){
            //find his neighbour
            for(std::size_t nei=0; nei<nodes_unsorted.size(); nei+2){
                
                //index of the element which is the neighbour (-1 if there is no neighbour)
                el_neighbour = -1;

                if( nei != el*N_nodes_per_element+n  && (nodes_unsorted[el*N_nodes_per_element+n]==nodes_unsorted[nei] 
                && nodes_unsorted[el*N_nodes_per_element+n+1]==nodes_unsorted[nei+1])
                ||(nodes_unsorted[el*N_nodes_per_element+n]==nodes_unsorted[nei+1] 
                && nodes_unsorted[el*N_nodes_per_element+n+1]==nodes_unsorted[nei])){

                    el_neighbour = (nei-(nei % (2*N_nodes_per_element))) / (2*N_nodes_per_element);  
                }
                /*   2 cases : - el_neighbour = -1 --> the values of u of the neighbour are in the limit conditions
                               - el_neighbour != -1 --> the values of u of the neighbour are in u 
                */
            }
            //calcul de la normale au bord
            gmsh::model::mesh::getNode(nodes_unsorted[el*N_nodes_per_element+n], nodeCoord1, nodeCoordParam1);
            gmsh::model::mesh::getNode(nodes_unsorted[el*N_nodes_per_element+n+1], nodeCoord2, nodeCoordParam2);

            // Computation of the normal. n = (-y , x)/(x^2+y^2)^(1/2)
            
            n_1 = nodeCoord1[1] - nodeCoord2[1]; // -y
            n_2 = nodeCoord2[0] - nodeCoord1[0]; // x

            double norm = sqrt(n_1 * n_1 + n_2 * n_2); // x^2 + y^2

            // Final norm.
            n_1 /= norm;
            n_2 /= norm;

            //sign of the normal (positive if the normal goes out of the element)
            gmsh::model::mesh::getNode(nodes_unsorted[el*N_nodes_per_element+((n+3)%N_nodes_per_element)], nodeCoord3, nodeCoordParam3);

            //if the dot product is negative then the normal has the good direction (sign --> positive)
            if(n_1*(nodeCoord3[0]-nodeCoord2[0]) + n_2*(nodeCoord3[1]-nodeCoord2[1])<0){
                
                sign=1;
            }else{
                sign=-1;
            }


            /*   2 cases : - el_neighbour = -1 --> the values of u of the neighbour are in the limit conditions
                        - el_neighbour != -1 --> the values of u of the neighbour are in u 
            */
            if(el_neighbour == -1){
                ///// besoin cond limites/////
                // 
                // u1 = boundary_cond(nodeCoord1, temps) // si jamais c'est dependant du temps 
                // u2 = boundary_cond(nodeCoord2, temps) 

            }else{
                //loop for the gauss points
                for(std::size_t gpt=0; gpt<N_gauss_point; gpt++){
                    //Loop for the L_j
                    for(std::size_t bf_j=0; bf_j < N_basis_function; bf_j++){
                        //loop for the L_i
                        for(std::size_t bf_i=0; bf_i < N_basis_function; bf_i++){

                            vectorF[el*N_basis_function + bf_i] += sign*(n_1*a_x + n_2*a_y) *integrationPoints[3+4*gpt]*bf[N_basis_function*gpt+ bf_i]*(1/2)*
                            (u[el*N_basis_function+bf_j] * bfg
                            - u[el_neighbour*N_basis_function+bf_j] * bfg);
                            //F[el*NumNodes + i] += sign*(n_1*a_x + n_2*a_y) * bf[NumNodes*g + i] * bf[NumNodes*g + j] * 0.5 * (u[e*NumNodes + j] - u[]
                    }
                }  
            }
        }
    }
}
