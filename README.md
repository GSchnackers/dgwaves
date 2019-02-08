# dgwaves

This projects shows how to use gmsh SDK on several platforms

## PC/windows/codeblocks

 * Go to the [gmsh web page](http://gmsh.info/) and download the latest version of gmsh-SDK: Windows-32bits (Code::Blocks provides a 32bits compiler)
 * Unzip it somewhere e.g. `C:\local\gmsh-4.1.4-Windows32-sdk`
 * Optional: you can create a soft link to shorten this path with the commands Open a terminal with `cmd.exe`, then:
```
c:
cd local 
mklink /J gmsh-sdk32 gmsh-4.1.4-Windows32-sdk
```
 * It is also very convenient to move the gmsh dll from the `lib` folder to the `bin` folder:
```
move gmsh-sdk32\lib\gmsh-4.1.dll gmsh-sdk32\bin
``` 
 * Overwrite `gmsh.h` with `gmsh.h_cwrap`:
```
copy /Y gmsh-sdk32\include\gmsh.h_cwrap gmsh-sdk32\include\gmsh.h
```
 * Add the following environment variables to your windows session:
```
set GMSHSDK=C:\local\gmsh-sdk32
set PATH=%GMSHSDK%\bin;%PATH%
set INCLUDE=%GMSHSDK%\include;%INCLUDE%
set LIB=%GMSHSDK%\lib;%LIB%
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%
``` 
In other words, you should add the `bin` folder of `gmsh-sdk32` to your `PATH` so that you can type `gmsh` in any terminal. The `bin` folder also contains the gmsh dynamic library, that must be found when `gmsh` (or your future solver) is run. The `INCLUDE` and `LIB` variables allows `CMake` to find the header `gmsh.h` and the library `gmsh.lib` respectively. Eventually, you may want to add the `lib` folder to your `PYTHONPATH` so that the command `import gmsh` works in python (it could be useful if you want to run python examples)
 * Download the latest version of CMake from the [CMake website](https://cmake.org/download/) and install it (add `cmake` to your `PATH` when asked, so that you can type `cmake` from a terminal)
 * Download [Code::Blocks](http://www.codeblocks.org/downloads). Choose the version with the mingw compiler included. Install it with the default options.
 * Add the MinGW compiler to your `PATH`
 * Clone this repository
 * ...


## PC/Linux/gcc

gmsh SDK requires libgfortran.so.3 to be installed on the system:
```
sudo apt-get install libgfortran3
```

## macOS/clang

