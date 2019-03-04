SET(SRCS main.cpp mysorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})