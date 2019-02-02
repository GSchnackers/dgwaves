@echo off
:: Compilation with msvc

set GMSHSDK=F:\local\gmsh-git-Windows64-sdk

set PATH=%GMSHSDK%\bin;%GMSHSDK%\lib;%PATH%
set INCLUDE=%GMSHSDK%\include;%INCLUDE%
set LIB=%GMSHSDK%\lib;%LIB%
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%

%comspec% /K ""C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" amd64"
