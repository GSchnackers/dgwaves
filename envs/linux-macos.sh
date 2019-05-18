# This file setup an environment which allows you to compile the code
# on Linux or macOS using the default compiler (gcc or clang)
#
# HOW TO USE THIS FILE?
#   open a terminal
#   . ./devenv-linux-macos.sh     # the first dot is important!
#   mkdir build
#   cd build
#   cmake ..
#   make
#   [executables are built in the bin/ folder]
#   ctest

# set the location of gmsh SDK ( **MODIFY THIS LINE FOR YOUR SYSTEM** )
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/media/nicolas/Seagate/Documents/Cours/Master1/Q2/integrated_project/dgwaves/gmsh-sdk/lib
export LD_LIBRARY_PATH
