@echo off
:: This file opens a terminal which allows you to compile the code with
:: a 64-bit mingw compiler
::     https://mingw-w64.org/
::     https://www.msys2.org/
::
:: HOW TO USE THIS FILE?
::   [check the PATHs below]
::   [run this file]
::   mkdir build
::   cd build
::   cmake -G "MinGW Makefiles" ..
::   mingw32-make
::   ctest
::
:: How to clean the "build" folder using cmd line?
::   cd build
::   rd /q /s .

echo setting MinGW64 environment...

:: set the location of the linear solver eigen.
set INCLUDE=C:\Users\Programme\eigen;%INCLUDE%

:: set the location of the header of the .h of the dgfem solver.
set INCLUDE=%CD%\Header;%INCLUDE%

:: set the location of the header of openmp.
set INCLUDE=C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\lib\gcc\x86_64-w64-mingw32\8.1.0\include;%INCLUDE%

:: set the location of the library of openmp.
set LIB=C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\lib\gcc\x86_64-w64-mingw32\8.1.0;%LIB%

:: set the location of gmsh SDK ( **MODIFY THIS LINE FOR YOUR SYSTEM** )
set GMSHSDK=C:\Users\Programme\gmsh-4.1.5-Windows64-sdk

:: where is gmsh.exe and gmsh-**.dll ? (HINT: copy gmsh-**.dll to the bin folder)
set PATH=%GMSHSDK%\bin;%PATH%

:: where is gmsh.h ? (rename gmsh.h_cwrap => gmsh.h)
set INCLUDE=%GMSHSDK%\include;%INCLUDE%

:: where is gmsh.lib ?
set LIB=%GMSHSDK%\lib;%LIB%

:: where is gmsh.py ? (required only if you want to use the python API)
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%

::set PATH=C:\mingw-w64\mingw64\bin;%PATH%
set PATH=C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin;%PATH%

:: open terminal
CD /d "%~dp0"
CD ..
%comspec% /K
