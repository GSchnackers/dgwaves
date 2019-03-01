SET(SRCS myedges.cpp mysorting.cpp my2Dnormal.cpp gaussIntegration.cpp)
ADD_EXECUTABLE(myedges ${SRCS})
TARGET_LINK_LIBRARIES(myedges ${GMSH_LIBRARIES})

SET(SRCS myview.cpp)
ADD_EXECUTABLE(myview ${SRCS})
TARGET_LINK_LIBRARIES(myview ${GMSH_LIBRARIES})
