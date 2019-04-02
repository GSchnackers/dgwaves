SET(SRCS main.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp invert.cpp Forward_Euler_method.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(main ${EIGEN_LIBRARIES})

SET(SRCS DG.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp initialCondition.cpp boundaryConditions.cpp invert.cpp Forward_Euler_method.cpp)
ADD_EXECUTABLE(DG ${SRCS})
TARGET_LINK_LIBRARIES(DG ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(DG ${EIGEN_LIBRARIES})

SET(SRCS mainTransport2D.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp initialCondition.cpp boundaryConditions.cpp invert.cpp Forward_Euler_method.cpp)
ADD_EXECUTABLE(mainTransport2D ${SRCS})
TARGET_LINK_LIBRARIES(mainTransport2D ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(mainTransport2D ${EIGEN_LIBRARIES})