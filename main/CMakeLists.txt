SET(SRCS main.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})

SET(SRCS DG.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp initialCondition.cpp)
ADD_EXECUTABLE(DG ${SRCS})
TARGET_LINK_LIBRARIES(DG ${GMSH_LIBRARIES})