#!/bin/sh

  module load cmake/3.11.1
  module load gcc/4.9.2

export CC=gcc
export CXX=g++

cd gmsh-sdk/
export FC=gfortran
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH}
cd ../

cd eigen/
export INCLUDE=${PWD}:${INCLUDE}

cd ../
rm -rf build/
mkdir build

cd build/
cmake ..
make
