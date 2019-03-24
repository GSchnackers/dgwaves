SET(SRCS main.cpp Initialization.cpp frontierCreation.cpp sortingNeighbouring.cpp getRealGradient.cpp invert.cpp getJacobianInverse.cpp normals.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})