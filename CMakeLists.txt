PROJECT(VTKLITE CXX Fortran)
# "Fortran" is required here because of MUMPS static library
cmake_minimum_required(VERSION 3.0)

# put all generated files in a single folder (required in windows for exe linked to dlls)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# build type is "" by default in Linux
IF(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# enable C++11
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" STD11CHECK)
IF(STD11CHECK)
    ADD_DEFINITIONS(-std=c++11 )
ENDIF()
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    ADD_DEFINITIONS(-Wno-literal-suffix)
ENDIF()

# find gmsh-sdk
# gmsh.h
find_path(GMSH_INCLUDE_DIRS NAMES "gmsh.h")
message(STATUS "GMSH_INCLUDE_DIRS=" ${GMSH_INCLUDE_DIRS})
if(NOT GMSH_INCLUDE_DIRS)
    message(FATAL_ERROR "gmsh.h not found!")
endif()
include_directories(${GMSH_INCLUDE_DIRS})

# libgmsh.so
find_library(GMSH_LIBRARIES gmsh)
message(STATUS "GMSH_LIBRARIES=" ${GMSH_LIBRARIES})
if(NOT GMSH_LIBRARIES)
    message(FATAL_ERROR "gmsh library not found!")
endif()


# program
SET(SRCS explore.cpp)
ADD_EXECUTABLE(explore ${SRCS})
target_link_libraries(explore ${GMSH_LIBRARIES})