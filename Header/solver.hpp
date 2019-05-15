// This header concers the functions related to the solver of the DG-FEM

#ifndef SOLVER_H
#define SOLVER_H
#include "structures.hpp"

/*
   Function that computes tha analytical solution
*/
 void compare(std::vector<double> & error, std::vector<double> & errorNodes, const Quantity & u,\
             const std::vector<double> & coordinates, const Element & mainElement,\
             const Simulation & simulation, const double mytime, const Properties & matProp);

/*
   Function that deals with the computation of the coefficient "k" useful for the Euler or Rugen-Kutta method.
*/
void computeCoeff(const Element & mainElement, const Element & frontierElement, const Simulation & simulation,\
                  const Properties & matProp, const double t, Quantity & u, Quantity & flux, \
                  std::vector<double> & k);

/*
   Function that initializes all quantities, parameters, etc... required for the simulation to work.
*/
void numericalInitializer(const Element & mainElement, Element & frontierElement, \
                          const Simulation & simulation, const PhysicalGroups & physicalGroups,\
                          Quantity & u, Quantity & flux, Properties & matProp);

/*
   Function that computes the Lax-Friederichs numerical flux for electromagnetism. "alpha" is a constant that
   weight the upwind and averaged part of the numerical flux.
*/
void numFluxELM(const Element & frontierElement, const Properties & matProp, const double alpha, Quantity & u, \
                Quantity & flux);

/*
   Function that computes the integration of the scalar product between the outside normals to elements and the 
   numerical fluxes over the frontier of the elements.
*/
void numFluxIntegration(const Quantity & flux, const Element & mainElement, const Element & frontierElement,\
                        std::vector<double> & fluxVector, int uNum);

/*
   Function that computes the simple upwind numerical flux in the case of a scalar transport.
*/
void numFluxUpwind(const Element & frontierElement, Quantity & flux);

/*
   Function that computes the physical flux cu at the nodes of the main elements and at the gauss points of the 
   frontier elements. "c" is the propagation speed of the quantity "u".
*/
void physFluxCu(const Quantity & u, const Element & mainElement, const Element & frontierElement,\
                Quantity & flux, std::vector<double> c);

/*
   Function tha computes the Lax-Friederichs flux for electromagnetism. 
*/
void physFluxELM(const Quantity & u, const Element & frontierElement, const Element & mainElement,\
                 const Properties & matProp, double t, Quantity & flux);

/*
   Function that set the boundary conditions at the nodes where they apply.
*/
void setBoundaryCondition(Element & frontierElement, const Simulation & simulation, \
                          const PhysicalGroups & physicalGroups, Quantity & u);

/*
   Function that computes the product between the stifness matrix and the nodal physical flux. The result is
   stored in "prod".
*/
void stiffnessFluxProd(const Element & mainElement, const Quantity & flux, std::vector<double> & prod, int uNum);

/*
   Function that handle the resolution of the DG-FEM method.
*/
void solver(const Element & mainElement, Element & frontierElement, const PhysicalGroups & physicalGroups,\
            View & view1, View & view2, Simulation & simulation);

/*
   Functions that computes the value of the quantity "u" at the gauss points. "force" forces the computation 
   at the outside of the domain.
*/
void valGp(Quantity & u, const Element & mainElement, const Element & frontierElement, int numU, \
           const Properties & matProp, double t = 0);

/*
   Write vector error in a file
*/
void writeError(const std::vector<double> & error, const Simulation & simulation);

#endif