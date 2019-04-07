SET(SRCS main.cpp Initialization.cpp frontierCreation.cpp sortingNeighbouring.cpp invert.cpp getJacobianInverse.cpp normals.cpp matrixMaker.cpp meshLoader.cpp boundaryCondition.cpp correspondance.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})