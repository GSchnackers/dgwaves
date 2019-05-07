// This header contains all function prototypes related to meshing.

#ifndef MESHING_H
#define MESHING_H

#include "structures.hpp"

/*
    Function that computes a vector containing the index correspondance between two adjacent nodes of two
    neighbours elements. It contains the local indices corresponding to the adjacent nodes.
*/
void correspondance(const Element & mainElement, Element & frontierElement);

/*
   Creates the elements making the frontier elements of the simulation. "mainElement" and "frontierElement" are
   respectively the  the elements of the bulk of the mesh and the elements between the elements in the bulk
   of the mesh. "meshDim" is the dimension of the bulk elements. "sortedNodes" is a vector of the nodes sorted of 
   the mesh, that is where each group of nodes making at the frontier of two main elements appear only once.
*/
void frontierCreation(const Element mainElement, Element & frontierElement, const int meshDim,\
                      const std::vector<int> sortedNodes);

/*
   Function that inverses the jacobian of the elements represented by "element".
*/
void getJacobiansInverse(Element & element);

/*
   Function that initializes an element. "element" is the element to be initialized. "gaussType" is the type of
   Gauss integration. "frontier" must be true if a frontier element is initialized. Else, it is false.
*/
void Initialization(Element & element, const int meshDim, std::string gaussType, bool frontier = false);

/* 
   Functions that invert a matrix. "matrix" contains the matrix to be inverted as a vector of the form
   [a11, a12, ..., a1N, a21, a22, ...] while "inverse" contains the inverse of that matrix also as a vector.
*/
void invert(std::vector<double> matrix, std::vector<double> & inverse);

/*
   This functions allow to compute the mass and stiffness matrices of the elements "element" required for
   the DG_FEM. "matrixType" is Sx, Sy, Sz, or M respectively for the stifness matrix with respect to x, y, or z 
   or the mass matrix.
*/
void matrixMaker(Element & element, std::string matrixType);

/* 
   Function that organizes the loading of the mesh of the simulation. "gaussType" is the GMSH integration Gauss 
   type. "physicalGroups" represent all information required on the physical group for the simulation to work.
*/
void meshLoader(Element & mainElements, Element & frontierElement, std::string & gaussType, \
                PhysicalGroups & physicalGroups, int mainDim);
/*
   Function that compute the normal to the frontier elements. The normal are always exterior to the first
   main element neighbour of the frontier element.
*/
void normals(Element & frontierElement, Element & mainElement);

/*
   Functions that produce the vector "nodeSorted" while determining the neighbours of each elements. 
*/
void sortingNeighbouring(const Element & mainElement, Element & frontierElement,\
                         std::vector<int> & nodeSorted);

#endif