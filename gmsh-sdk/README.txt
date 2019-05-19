This is the binary Software Development Kit (SDK) for Gmsh 4.1.5:

  * Operating system: Linux64-sdk (Linux)
  * C++ compiler: /usr/bin/c++
  * C++ compiler ID: GNU
  * C++ compiler version: 4.7
  * C++ compiler flags: -std=c++11 -O2 -g -DNDEBUG
  * Build options: 64Bit Ann Bamg Bfgs Blas Blossom Cgns DIntegration Dlopen DomHex Fltk Gmm Hxt Hxt3D Jpeg[fltk] Kbipack Lapack LinuxJoystick MathEx Med Mesh Metis Mmg3d Mpeg Netgen ONELAB ONELABMetamodel OpenCASCADE OpenGL OptHom PETSc Parser Plugins Png[fltk] Post QuadTri SLEPc Solver TetGen/BR Voro++ Zlib

Gmsh is distributed under the terms of the GNU General Public License: see
share/doc/gmsh/LICENSE.txt and share/doc/gmsh/CREDITS.txt. For additional Gmsh
resources, see http://gmsh.info.

SDK layout:

  * lib/*gmsh*.{so,dll,dylib}: shared Gmsh library
  * lib/gmsh.lib: import library (Windows only)
  * lib/gmsh.py: Python module
  * lib/gmsh.jl: Julia module
  * include/gmsh.h: C++ API header
  * include/gmshc.h: C API header
  * include/gmsh.h_cwrap: C++ wrapper of the C API (see the `Notes' below)
  * bin/gmsh: gmsh executable (linked with the shared Gmsh library)
  * share/doc/gmsh/demos/api : API examples in C++, C, Python and Julia

Notes:

  * The C API should work with most compilers.

  * The C++ API will only work if your compiler has the same Application Binary
    Interface (ABI) as the compiler used to build this SDK.

    For example, the Linux SDK is currently compiled with GCC 4. You would need
    to specify "-D_GLIBCXX_USE_CXX11_ABI=0" if you use GCC 5 or newer.

  * If your C++ compiler does not provide the same ABI as the C++ compiler used
    to build the SDK and if there are no compatibility flags available, you can
    rename `gmsh.h_cwrap' as `gmsh.h': this implementation redefines the C++ API
    in terms of the C API. Using this header will lead to (slightly) reduced
    performance compared to using the native Gmsh C++ API from the original
    `gmsh.h', as it entails additional data copies between this C++ wrapper, the
    C API and the native C++ code.

    For example, the Windows SDK is currently compiled using GCC 5. To compile a
    C++ example with Microsoft Visual Studio 2017 in the Visual Studio shell and
    run it, you would do:

    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h_cwrap gmsh.h
    C:\gmsh-git-Windows64-sdk> cl /Iinclude share\doc\gmsh\demos\api\simple.cpp lib\gmsh.lib
    C:\gmsh-git-Windows64-sdk> cd lib
    C:\gmsh-git-Windows64-sdk\lib> ..\simple.exe

  * The shared Gmsh library references several other shared system libraries. On
    Linux for example, it depends on libgfortran.so.3: you will need to install
    this library (search your package manager for it) if your Linux distribution
    does not provide it by default, or if it comes with a different version.
