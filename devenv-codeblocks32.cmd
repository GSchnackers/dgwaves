@echo off
:: Compilation with Code::Blocks (+ embedded mingw (32bits))
::
:: cmake -G "CodeBlocks - MinGW Makefiles" ..

set GMSHSDK=C:\local\gmsh-4.1.4-Windows32-sdk

set PATH=%GMSHSDK%\bin;%GMSHSDK%\lib;%PATH%
set INCLUDE=%GMSHSDK%\include;%INCLUDE%
set LIB=%GMSHSDK%\lib;%LIB%
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%

%comspec% /K "C:\Program Files (x86)\CodeBlocks\MinGW\mingwvars.bat"
