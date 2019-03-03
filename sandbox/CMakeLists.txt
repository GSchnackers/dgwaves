SET(SRCS myedges.cpp mysorting.cpp normal.cpp gaussIntegration.cpp gradient.cpp)
ADD_EXECUTABLE(myedges ${SRCS})
TARGET_LINK_LIBRARIES(myedges ${GMSH_LIBRARIES})

SET(SRCS myview.cpp)
ADD_EXECUTABLE(myview ${SRCS})
TARGET_LINK_LIBRARIES(myview ${GMSH_LIBRARIES})
