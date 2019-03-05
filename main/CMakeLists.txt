<<<<<<< HEAD
SET(SRCS main.cpp mysorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp)
=======
SET(SRCS main.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp)
>>>>>>> 28274124b591580bdd9310cf9a19343e1a3dc581
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})

SET(SRCS DG.cpp sorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp neighbours.cpp)
ADD_EXECUTABLE(DG ${SRCS})
TARGET_LINK_LIBRARIES(DG ${GMSH_LIBRARIES})