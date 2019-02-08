# dgwaves

This projects shows how to use gmsh SDK on several platforms

## PC/windows/codeblocks

 * Go to the [gmsh web page](http://gmsh.info/) and download the latest version of gmsh-SDK: Windows-64bits (or Windows-32bits if you own an very old PC)
 * Unzip it somewhere e.g. `C:\local\gmsh-4.1.4-Windows64-sdk`
 * Optional: you can create a soft link to shorten this path with the commands Open a terminal with `cmd.exe`, then:
```
c:
cd local 
mklink /J gmsh-sdk gmsh-4.1.4-Windows64-sdk
```
 * It is also convenient to move the gmsh dll from the `lib` folder to the `bin` folder:
```
move gmsh-sdk\lib\gmsh-4.1.dll gmsh-sdk\bin
``` 
 * Add the following environment variables to your windows session:
```
set GMSHSDK=C:\local\gmsh-sdk
set PATH=%GMSHSDK%\bin;%PATH%
set INCLUDE=%GMSHSDK%\include;%INCLUDE%
set LIB=%GMSHSDK%\lib;%LIB%
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%
``` 
In other words, you should add the `bin` and `lib` folder of gmsh-sdk to your `PATH` so that you can type `gmsh` in any terminal. The `lib` folder contains the gmsh dynamic library, that must be found when `gmsh` (or your future solver) is run. The `INCLUDE` and `LIB` variables allows `CMake` to find the header `gmsh.h` and the library `gmsh.lib`. Eventually, you may want to add the `lib` folder to your `PYTHONPATH` so that the command `import gmsh` works in python (it could be useful if you want to run python examples)
 * Download the latest version of CMake from the [CMake website](https://cmake.org/download/) and install it.
 *


## PC/Linux/gcc

gmsh SDK requires libgfortran.so.3 to be installed on the system:
```
sudo apt-get install libgfortran3
```

## macOS/clang

