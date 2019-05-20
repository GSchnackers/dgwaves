#!/bin/sh

module load cmake/3.11.1
module load gcc/4.9.2

export LD_LIBRARY_PATH=gmsh-sdk/lib
export FC=gfortran

rm -rf build/  
mkdir build

cd build/
cmake .. 
make

cd ..
