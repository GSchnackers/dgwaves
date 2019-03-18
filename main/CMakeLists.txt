SET(SRCS main.cpp Initialization.cpp frontierCreation.cpp sortingNeighbouring.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})