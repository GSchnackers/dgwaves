<<<<<<< HEAD
SET(SRCS main.cpp mysorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp)
=======
SET(SRCS main.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp)
>>>>>>> 0db2644af05e198df9e033cd3ae33365b9c64cf3
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})