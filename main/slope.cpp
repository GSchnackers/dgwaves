
#include <cstdio>
#include <iostream>
#include <vector>
#include <gmsh.h>
#include "functions.h"

void slope(std::vector<int> & nodeTags2D, std::vector<int> & nodeTags2DPlusBC, std::vector<double> & nodeCoord,\
            std::vector<double> & nodeCoordParam, double mytime, double value, std::vector<double> & matrixF,\
            std::vector<double> & uPlusBC, std::vector<double> & u, std::vector<double> & matrixS,\
            std::vector<int> & tagElement1DSorted, std::vector<int> & upwind, std::vector<int> & neighbours1D,\
            std::vector<int> & indicesNei1, std::vector<int> & indicesNei2, std::vector<double> & matrixM_Inverted,\
            std::vector<double> & dudt, std::vector<int> & elementTags2D, int numNodes2D, int NumNodesSide){

        // BC
        for(std::size_t i=nodeTags2D.size(); i<nodeTags2DPlusBC.size(); i++){
            gmsh::model::mesh::getNode(nodeTags2DPlusBC[i], nodeCoord, nodeCoordParam);
            boundaryConditions(nodeCoord, mytime, value);
            uPlusBC[i] = value;
        }

        // TEST
        /*
        for(size_t i = 0; i < uPlusBC.size(); i++){
            std::cout << "uPlusBC[" << std::to_string(i) << "]  : " << std::to_string(uPlusBC[i]) << " at mytime " << std::to_string(mytime) << "\n";
        }
        */

        
        // declaration vector F (time dependent)
        std::vector<double> vectorF(nodeTags2D.size(), 0.0);

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
                    
                }
                if(neighbours1D[ed*2 + 1] != -1){
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
                    
                }
                if(neighbours1D[ed*2 + 1] != -1){
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


        // TEST vector F at each time
        /*
        std::cout <<"time = " << std::to_string(mytime) << "\n";
        
        for(std::size_t el=0; el < elementTags2D.size(); el++){
            std::cout << "el " << std::to_string(el) << "\n";
            //loop for i of F_{ij}
            for(std::size_t i=0; i < numNodes2D; i++){
                
                std::cout << std::to_string(vectorF[el*numNodes2D + i]) << " ";
                
            std::cout << "\n";
            }
        std::cout << "\n";
        }
        */

        // declaration vector Su
        std::vector<double> VectorSu(elementTags2D.size()*numNodes2D, 0.0);

        //Computation of VectorSu = S.u
        for(std::size_t el=0; el<elementTags2D.size(); el++){
            for(std::size_t i=0; i<numNodes2D; i++){
                for(std::size_t j=0; j<numNodes2D; j++){

                    VectorSu[el*numNodes2D + i] += matrixS[el*numNodes2D*numNodes2D + i*numNodes2D + j] * u[el*numNodes2D + j];
                }
            }
        }

        // TEST Su
        /*
        for(std::size_t el=0; el < elementTags2D.size(); el++){
            std::cout << "el " << std::to_string(el) << "\n";
            //loop for i of F_{ij}
            for(std::size_t i=0; i < numNodes2D; i++){
                
                std::cout << "Su["<< std::to_string(el*numNodes2D + i) << "] = " << std::to_string(VectorSu[el*numNodes2D + i]) << "\n";
                
            }
        std::cout << "\n";
        }
        */

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
                                                    (VectorSu[el*numNodes2D + j] + vectorF[el*numNodes2D + j]);
                }
            }
        }
}